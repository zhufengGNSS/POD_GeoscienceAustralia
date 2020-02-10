MODULE m_pd_empirical


! ----------------------------------------------------------------------
! MODULE: m_pd_empirical
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the m_pd_empirical subroutine
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	24 September 2018
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE pd_empirical (rsat, vsat, GMearth, Yangle, frame, Femp, PDr, PDv, PD_param)


! ----------------------------------------------------------------------
! SUBROUTINE: pd_empirical.f03
! ----------------------------------------------------------------------
! Purpose:
!  Acceleration vector and Partial derivatives of the emprirical force model 
! ----------------------------------------------------------------------
! Input arguments:
! - rsat:			Satellite Position vector (m)   in inertial frame (ICRF)
! - vsat:			Satellite Velocity vector (m/s) in inertial frame (ICRF)
! 
! Output arguments:
! - Femp:			Acceleration vector cartesian components in inertial frame (ICRF)
! - PDr: 			Partial derivatives matrix of the acceleration w.r.t. the position vector in ICRF
! - PDv: 			Partial derivatives matrix of the acceleration w.r.t. the velocity vector in ICRF
! - PD_param: 		Partial derivatives matrix of the acceleration w.r.t. the (force-related) unknown parameters in ICRF
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
! 			Cooperative Research Centre for Spatial Information, Australia
! Created:	September 2018
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_d), INTENT(IN), DIMENSION(3) :: rsat, vsat
      REAL (KIND = prec_q), INTENT(IN) :: GMearth
      REAL (KIND = prec_q), INTENT(IN) :: Yangle
      INTEGER (KIND = prec_int2), INTENT(IN) :: frame
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_d), INTENT(OUT) :: Femp(3)
      REAL (KIND = prec_d), INTENT(OUT) :: PDr(3,3), PDv(3,3)
      REAL (KIND = prec_d), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: PD_param
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(3) :: rsat_icrf, vsat_icrf
      REAL (KIND = prec_d) :: fx, fy, fz
      REAL (KIND = prec_q) :: ax, ay, az
      INTEGER (KIND = prec_int8) :: i , j 
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: N_param, PD_Param_ID, Param_Bias, Param_CPR
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: PD_Bias_param, PD_CPR_param 
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(3) :: pd_aRTN_biasR, pd_aEMP_biasR 
      REAL (KIND = prec_d), DIMENSION(3) :: pd_aRTN_biasT, pd_aEMP_biasT
      REAL (KIND = prec_d), DIMENSION(3) :: pd_aRTN_biasN, pd_aEMP_biasN 
      REAL (KIND = prec_q) :: Bias_RTN(3), Bias_icrf(3)    
      REAL (KIND = prec_q) :: CPR_RTN(3), CPR_icrf(3)    
	  INTEGER (KIND = prec_int2) :: nCPR
      REAL (KIND = prec_q) :: CPR_Cr, CPR_Sr, CPR_Ct, CPR_St, CPR_Cn, CPR_Sn  
      REAL (KIND = prec_q) :: aCPR_r, aCPR_t, aCPR_n 
	  REAL (KIND = prec_q), DIMENSION(3) :: pd_aRTN_Cr, pd_aRTN_Sr, pd_aRTN_Ct, pd_aRTN_St, pd_aRTN_Cn, pd_aRTN_Sn  
      REAL (KIND = prec_q) :: pd_aR_u, pd_aT_u, pd_aN_u 
      REAL (KIND = prec_q) :: Rrtn(3,3), Rrtn_inv(3,3)
      REAL (KIND = prec_q) :: kepler(9), u_deg, u_rad
	  REAL (KIND = prec_q) :: dxyz, u_deg_dx, u_deg_dy, u_deg_dz, u_rad_dx, u_rad_dy, u_rad_dz
	  REAL (KIND = prec_q) :: Pd_aCPR_xyz(3,3), Pd_aCPR_rtn_xyz(3,3), PD_u_xyz(3), PD_u_x, PD_u_y, PD_u_z 
	  REAL (KIND = prec_q), DIMENSION(3) :: r_dx, r_dy, r_dz
	  REAL (KIND = prec_q), DIMENSION(9) :: kepler_dx, kepler_dy, kepler_dz
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: Rcrf_bff(3,3), Rrtn_bff(3,3)
! ----------------------------------------------------------------------

	  
Param_Bias = 0
Param_CPR  = 0

! ----------------------------------------------------------------------
! Partial derivatives w.r.t. unknown parameters
N_param = NPARAM_glb
If (N_param == 0) Then
N_param = 1
End If
ALLOCATE (PD_param(3,N_param), STAT = AllocateStatus)
! init needed
PD_param = 0.d0

! Bias partial derivatives matrix allocation
PD_Param_ID = 0
If (EMP_Bias_glb(1) == 1) Then
	PD_Param_ID = PD_Param_ID + 1
End IF
If (EMP_Bias_glb(2) == 1) Then
	PD_Param_ID = PD_Param_ID + 1
End IF
If (EMP_Bias_glb(3) == 1) Then
	PD_Param_ID = PD_Param_ID + 1
End IF
ALLOCATE (PD_Bias_param(3,PD_Param_ID), STAT = AllocateStatus)
Param_Bias = PD_Param_ID
PD_Bias_param = 0.d0
PD_Param_ID = 0

! CPR partial derivatives matrix allocation
PD_Param_ID = 0
If (EMP_CPR_glb(1) == 1) THEN		
	PD_Param_ID = PD_Param_ID + 2
End IF
If (EMP_CPR_glb(2) == 1) THEN		
	PD_Param_ID = PD_Param_ID + 2
End IF
If (EMP_CPR_glb(3) == 1) THEN		
	PD_Param_ID = PD_Param_ID + 2
End If
ALLOCATE (PD_CPR_param(3,PD_Param_ID), STAT = AllocateStatus)
Param_CPR = PD_Param_ID
PD_CPR_param = 0.d0
PD_Param_ID = 0

! ----------------------------------------------------------------------
! State Vector in ICRF
rsat_icrf = rsat
vsat_icrf = vsat
! ----------------------------------------------------------------------
 
! ----------------------------------------------------------------------
! Reference Frame of empirical forces:
! 1. Orbital Frame
! 2. Body-fixed Frame
! ----------------------------------------------------------------------
if (frame == 1) then

! Inertial (GCRF) to Orbital Frame
Call orb_frame(rsat_icrf, vsat_icrf, Rrtn)	

elseif (frame == 2) then

! Inertial (GCRF) to Body-fixed Frame
CALL crf_bff (rsat_icrf, vsat_icrf, Yangle, Rcrf_bff, Rrtn_bff)
Rrtn = Rcrf_bff

end if

! Orbital frame to Inertial frame (GCRF) : Inverse matrix
CALL matrix_inv3 (Rrtn, Rrtn_inv)
! ----------------------------------------------------------------------
	
! ----------------------------------------------------------------------
! Bias accelerations	  
If (EMP_Bias_glb(1) == 1) Then
	Bias_RTN(1) = Bias_accel_glb(1)
End IF
If (EMP_Bias_glb(2) == 1) Then
	Bias_RTN(2) = Bias_accel_glb(2)
End IF
If (EMP_Bias_glb(3) == 1) Then
	Bias_RTN(3) = Bias_accel_glb(3)
End IF
If (EMP_Bias_glb(1) == 1 .OR. EMP_Bias_glb(2) == 1 .OR. EMP_Bias_glb(3) == 1) Then
	!Bias_RTN = Bias_accel_glb
	Bias_icrf = MATMUL(Rrtn_inv, Bias_RTN)
Else
	Bias_icrf = (/ 0.0D0, 0.0D0, 0.0D0 /)
End IF
! ----------------------------------------------------------------------
! Partial Derivatives of Bias accelerations w.r.t. state vector:
! PD_Bias_r = 0
! PD_Bias_v = 0
! ----------------------------------------------------------------------
! Partial Derivatives of Bias accelerations w.r.t. unkown parameters i.e. Biases
! PD_Bias_param = [pd_aEMP_biasR pd_aEMP_biasT pd_aEMP_biasN] 
PD_Param_ID = 0
If (EMP_Bias_glb(1) == 1) Then
	pd_aRTN_biasR = (/ 1.0D0, 0.0D0, 0.0D0 /)
	pd_aEMP_biasR = MATMUL(Rrtn_inv, pd_aRTN_biasR)
	PD_Param_ID = PD_Param_ID + 1
	PD_Bias_param(1,PD_Param_ID) = pd_aEMP_biasR(1)
	PD_Bias_param(2,PD_Param_ID) = pd_aEMP_biasR(2)
	PD_Bias_param(3,PD_Param_ID) = pd_aEMP_biasR(3)
End IF
If (EMP_Bias_glb(2) == 1) Then
	pd_aRTN_biasT = (/ 0.0D0, 1.0D0, 0.0D0 /)
	pd_aEMP_biasT = MATMUL(Rrtn_inv, pd_aRTN_biasT)
	PD_Param_ID = PD_Param_ID + 1
	PD_Bias_param(1,PD_Param_ID) = pd_aEMP_biasT(1)
	PD_Bias_param(2,PD_Param_ID) = pd_aEMP_biasT(2)
	PD_Bias_param(3,PD_Param_ID) = pd_aEMP_biasT(3)
End IF
If (EMP_Bias_glb(3) == 1) Then
	pd_aRTN_biasN = (/ 0.0D0, 0.0D0, 1.0D0 /)
	pd_aEMP_biasN = MATMUL(Rrtn_inv, pd_aRTN_biasN)
	PD_Param_ID = PD_Param_ID + 1
	PD_Bias_param(1,PD_Param_ID) = pd_aEMP_biasN(1)
	PD_Bias_param(2,PD_Param_ID) = pd_aEMP_biasN(2)
	PD_Bias_param(3,PD_Param_ID) = pd_aEMP_biasN(3)
End IF
PD_Param_ID = 0
! ----------------------------------------------------------------------
!print *,"PD_Param_ID", PD_Param_ID
!print *,"PD_Bias_param(:,1)", PD_Bias_param(:,1)
!print *,"PD_Bias_param(:,2)", PD_Bias_param(:,2)
!print *,"PD_Bias_param(:,3)", PD_Bias_param(:,3)
!print *,"PD_Bias_param", PD_Bias_param
!print *,"Rrtn_inv", Rrtn_inv


! ----------------------------------------------------------------------
! Cycle-per-revolution accelerations	  
! ----------------------------------------------------------------------
If (EMP_CPR_glb(1) == 1 .OR. EMP_CPR_glb(2) == 1 .OR. EMP_CPR_glb(3) == 1) Then
	! CPR coefficients read through module mdl_param.f03
	nCPR = EMP_nCPR_glb
	CPR_Cr = CPR_CS_glb(1,1)
	CPR_Sr = CPR_CS_glb(1,2)
	CPR_Ct = CPR_CS_glb(2,1)
	CPR_St = CPR_CS_glb(2,2)
	CPR_Cn = CPR_CS_glb(3,1)
	CPR_Sn = CPR_CS_glb(3,2)
	! Keplerian element u: argument of latitude
	Call kepler_z2k (rsat_icrf, vsat_icrf, GMearth, kepler)
	u_deg = kepler(9)
	u_rad = u_deg * (PI_global / 180.D0)
! ----------------------------------------------------------------------
! Acceleration vector	
	If (EMP_CPR_glb(1) == 1) THEN		
		! Radial direction: acceleration vector
		aCPR_r = CPR_Cr * cos(nCPR*u_rad) + CPR_Sr * sin(nCPR*u_rad)
		! Partial derivative w.r.t u argument 
		pd_aR_u = -1.0D0 * CPR_Cr * sin(nCPR*u_rad) * nCPR + CPR_Sr * cos(nCPR*u_rad) * nCPR			! 777 777 777
		!pd_aR_u =          CPR_Cr * sin(nCPR*u_rad) * nCPR + CPR_Sr * cos(nCPR*u_rad) * nCPR
	else
		aCPR_r  = 0.0D0
		pd_aR_u	= 0.0D0	
	End If
	! Along-track direction
	If (EMP_CPR_glb(2) == 1) THEN		
		aCPR_t = CPR_Ct * cos(nCPR*u_rad) + CPR_St * sin(nCPR*u_rad)
		pd_aT_u = -1.0D0 * CPR_Ct * sin(nCPR*u_rad) * nCPR + CPR_St * cos(nCPR*u_rad) * nCPR	
	else
		aCPR_t  = 0.0D0	
		pd_aT_u	= 0.0D0	
	End If
	! Cross-track direction
	If (EMP_CPR_glb(3) == 1) THEN		
		aCPR_n = CPR_Cn * cos(nCPR*u_rad) + CPR_Sn * sin(nCPR*u_rad)
		pd_aN_u = -1.0D0 * CPR_Cn * sin(nCPR*u_rad) * nCPR + CPR_Sn * cos(nCPR*u_rad) * nCPR
	else
		aCPR_n  = 0.0D0
		pd_aN_u	= 0.0D0			
	End If
! Acceleration vector in orbital frame
	CPR_RTN = (/ aCPR_r, aCPR_t, aCPR_n /)
! Acceleration vector in inertial frame
	CPR_icrf = MATMUL(Rrtn_inv, CPR_RTN)
Else
	CPR_icrf = (/ 0.0D0, 0.0D0, 0.0D0 /)
! ----------------------------------------------------------------------
! dummy init until more knowledge known
        u_rad = 0.d0
        pd_aT_u = 0.d0
        pd_aN_u = 0.d0
        pd_aR_u = 0.d0
        nCPR = 0
End IF	
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Partial derivatives w.r.t state vector
! ----------------------------------------------------------------------
! PDr, PDv
! PDr = [PD_aCPRx_x PD_aCPRx_y PD_aCPRx_z 
! 		 PD_aCPRy_x	PD_aCPRy_y PD_aCPRy_z
!		 PD_aCPRz_x PD_aCPRz_y PD_aCPRz_z] 
! ----------------------------------------------------------------------
! d(u)/d(xyz) derivatives
dxyz = 0.001D0
r_dx = (/rsat_icrf(1)+dxyz, rsat_icrf(2), rsat_icrf(3) /)
r_dy = (/rsat_icrf(1), rsat_icrf(2)+dxyz, rsat_icrf(3) /)
r_dz = (/rsat_icrf(1), rsat_icrf(2), rsat_icrf(3)+dxyz /)

Call kepler_z2k (r_dx, vsat_icrf, GMearth, kepler_dx)
Call kepler_z2k (r_dy, vsat_icrf, GMearth, kepler_dy)
Call kepler_z2k (r_dz, vsat_icrf, GMearth, kepler_dz)

u_deg_dx = kepler_dx(9)
u_rad_dx = u_deg_dx * (PI_global / 180.D0)
u_deg_dy = kepler_dy(9)
u_rad_dy = u_deg_dy * (PI_global / 180.D0)
u_deg_dz = kepler_dz(9)
u_rad_dz = u_deg_dz * (PI_global / 180.D0)

PD_u_x = (u_rad_dx - u_rad) / dxyz
PD_u_y = (u_rad_dy - u_rad) / dxyz
PD_u_z = (u_rad_dz - u_rad) / dxyz
PD_u_xyz = (/ PD_u_x, PD_u_y, PD_u_z /)
! ----------------------------------------------------------------------
! Pd_aCPR_rtn_xyz
!Pd_aCPR_rtn_xyz (1,1:3) = MATMUL(pd_aR_u, PD_u_xyz)
Pd_aCPR_rtn_xyz (1,1) = pd_aR_u * PD_u_x
Pd_aCPR_rtn_xyz (1,2) = pd_aR_u * PD_u_y
Pd_aCPR_rtn_xyz (1,3) = pd_aR_u * PD_u_z

!Pd_aCPR_rtn_xyz (2,1:3) = MATMUL(pd_aT_u, PD_u_xyz)
Pd_aCPR_rtn_xyz (2,1) = pd_aT_u * PD_u_x
Pd_aCPR_rtn_xyz (2,2) = pd_aT_u * PD_u_y
Pd_aCPR_rtn_xyz (2,3) = pd_aT_u * PD_u_z

!Pd_aCPR_rtn_xyz (3,1:3) = MATMUL(pd_aN_u, PD_u_xyz)
Pd_aCPR_rtn_xyz (3,1) = pd_aN_u * PD_u_x
Pd_aCPR_rtn_xyz (3,2) = pd_aN_u * PD_u_y
Pd_aCPR_rtn_xyz (3,3) = pd_aN_u * PD_u_z

! Transformation to ICRF
Pd_aCPR_xyz = MATMUL(Rrtn_inv, Pd_aCPR_rtn_xyz)
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! Partial derivatives w.r.t unknown parameters to be estimated
PD_Param_ID = 0
If (EMP_CPR_glb(1) == 1) THEN		
! Radial component
! C term
	pd_aRTN_Cr = (/ cos(nCPR*u_rad), 0.0D0, 0.0D0 /) 
	PD_Param_ID = PD_Param_ID + 1
	PD_CPR_param(1:3,PD_Param_ID) = MATMUL(Rrtn_inv, pd_aRTN_Cr)
! S term
	pd_aRTN_Sr = (/ sin(nCPR*u_rad), 0.0D0, 0.0D0 /)
	PD_Param_ID = PD_Param_ID + 1
	PD_CPR_param(1:3,PD_Param_ID) = MATMUL(Rrtn_inv, pd_aRTN_Sr)
	
End IF
If (EMP_CPR_glb(2) == 1) THEN		
! Along-track component
! C term
	pd_aRTN_Ct = (/ 0.0D0, cos(nCPR*u_rad), 0.0D0 /) 
	PD_Param_ID = PD_Param_ID + 1
	PD_CPR_param(1:3,PD_Param_ID) = MATMUL(Rrtn_inv, pd_aRTN_Ct)	
! S term
	pd_aRTN_St = (/ 0.0D0, sin(nCPR*u_rad), 0.0D0 /)
	PD_Param_ID = PD_Param_ID + 1
	PD_CPR_param(1:3,PD_Param_ID) = MATMUL(Rrtn_inv, pd_aRTN_St)	
End IF
If (EMP_CPR_glb(3) == 1) THEN		
! Cross-track component
! C term	
	pd_aRTN_Cn = (/ 0.0D0, 0.0D0, cos(nCPR*u_rad) /) 
	PD_Param_ID = PD_Param_ID + 1
	PD_CPR_param(1:3,PD_Param_ID) = MATMUL(Rrtn_inv, pd_aRTN_Cn)	
! S term	
	pd_aRTN_Sn = (/ 0.0D0, 0.0D0, sin(nCPR*u_rad) /)
	PD_Param_ID = PD_Param_ID + 1
	PD_CPR_param(1:3,PD_Param_ID) = MATMUL(Rrtn_inv, pd_aRTN_Sn)		
End If	
! ----------------------------------------------------------------------

! End of CPR terms
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Summary of Empirical accelerations	
Femp = Bias_icrf + CPR_icrf
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Partial derivatives w.r.t state vector
PDr = Pd_aCPR_xyz
!PDv
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Partial derivatives w.r.t unknown parameters to be estimated
! FIXME: tried but still cannot get rid of warnings in this block
! ----------------------------------------------------------------------
!PD_param = [PD_Bias_param PD_CPR_param]
If (Param_Bias>0 .and. Param_CPR==0) Then
	PD_param = PD_Bias_param
Else If (Param_Bias==0 .and. Param_CPR>0) Then
	PD_param = PD_CPR_param
Else If (Param_Bias>0 .and. Param_CPR>0) Then
i = 0
j = 0
Do i = 1 , 3
Do j = 1 , Param_Bias
	PD_param(i,j) = PD_Bias_param(i,j)
End Do
End Do

i = 0
j = 0
Do j = 1 , Param_CPR
Do i = 1 , 3
	PD_param(i,j+Param_Bias) = PD_CPR_param(i,j)
End Do
End Do
End IF
! ----------------------------------------------------------------------
i = size(PD_param, DIM = 1)
j = size(PD_param, DIM = 2)
!print *,"PD_param DIM",i,j
!print *,"Param_Bias,Param_CPR",Param_Bias,Param_CPR
!print *,"PD_param",PD_param

! ----------------------------------------------------------------------
! FIXME: calculate PDv. Initialise to zero for now
! ----------------------------------------------------------------------
PDv = 0.d0

END SUBROUTINE

END Module
