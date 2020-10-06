MODULE m_pd_pulses


! ----------------------------------------------------------------------
! MODULE: m_pd_pulses
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the pd_pulses subroutine
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou at Geoscience Australia
! Created:	17 August 2020
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE pd_pulses (rsat, vsat, mjd_t, delta_v, mjd_ti, dir, Fpulse, PDr, PDv, PD_param)


! ----------------------------------------------------------------------
! SUBROUTINE: pd_empirical.f03
! ----------------------------------------------------------------------
! Purpose:
!  Acceleration vector and Partial derivatives of pseudo-stochastic pulses (Velocity changes) 
! ----------------------------------------------------------------------
! Input arguments:
! - rsat:			Satellite Position vector (m)   in inertial frame (ICRF)
! - vsat:			Satellite Velocity vector (m/s) in inertial frame (ICRF)
! 
! Output arguments:
! - Fpulse:		Acceleration vector cartesian components in inertial frame (ICRF)
! - PDr: 			Partial derivatives matrix of the acceleration w.r.t. the position vector in ICRF
! - PDv: 			Partial derivatives matrix of the acceleration w.r.t. the velocity vector in ICRF
! - PD_param: 		Partial derivatives matrix of the acceleration w.r.t. the (force-related) unknown parameters in ICRF
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou at Geoscience Australia
! Created:	17 August 2020
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
      INTEGER (KIND = prec_int2), INTENT(IN) :: dir
      REAL (KIND = prec_d), INTENT(IN) :: delta_v(3)
      REAL (KIND = prec_d), INTENT(IN) :: mjd_t, mjd_ti
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_d), INTENT(OUT) :: Fpulse (3)
      REAL (KIND = prec_d), INTENT(OUT) :: PDr(3,3), PDv(3,3)
      REAL (KIND = prec_d), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: PD_param
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(3) :: rsat_icrf, vsat_icrf
      REAL (KIND = prec_q) :: er(3), et(3), en(3), e_unit_dir(3) 
      REAL (KIND = prec_q) :: Rrtn(3,3), Rrtn_inv(3,3)
      REAL (KIND = prec_d) :: delta_dirac
      REAL (KIND = prec_d) :: delta_v_ti
      INTEGER (KIND = prec_int2) :: N_param_pulses
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
!      REAL (KIND = prec_q) :: Rrtn(3,3), Rrtn_inv(3,3)
!      REAL (KIND = prec_d) :: Rcrf_bff(3,3), Rrtn_bff(3,3)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
!N_param_pulses = N_PULSE_param_glb
!ALLOCATE (PD_param(3,N_param_pulses), STAT = AllocateStatus)
ALLOCATE (PD_param(3,1), STAT = AllocateStatus)
PD_param = 0.0D0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! State Vector in ICRF
rsat_icrf = rsat
vsat_icrf = vsat
! ----------------------------------------------------------------------
 
! ----------------------------------------------------------------------
! Unit vecotr of the direction of the velocity pulse:
! ----------------------------------------------------------------------
! 1. Radial direction 
! 2. Along-track
! 3. Cross-track
! ----------------------------------------------------------------------
CALL orb_frame2_unit(rsat_icrf, vsat_icrf, Rrtn, er, et, en)

! Orbital frame to Inertial frame (GCRF) : Inverse matrix
CALL matrix_inv3 (Rrtn, Rrtn_inv)

if (dir == 1) then
! Radial direction
e_unit_dir = er
delta_v_ti = delta_v(1)

elseif (dir == 2) then
! Tangential direction
e_unit_dir = et
delta_v_ti = delta_v(2)

elseif (dir == 3) then
! Normal direction
e_unit_dir = en
delta_v_ti = delta_v(3)

end if
! ----------------------------------------------------------------------
	

! ----------------------------------------------------------------------
! Dirac's "delta" function at the current epoch t (mjd_t)
! ----------------------------------------------------------------------
delta_dirac = 0.0D0

!IF (t == ti) THEN
IF ( abs(mjd_t - mjd_ti) < 1.0D-06 ) THEN
	delta_dirac = 1.0D0
END IF
! ----------------------------------------------------------------------

	
! ----------------------------------------------------------------------
! Summary of Pseudo-stochastic pulses	
! ----------------------------------------------------------------------
!Fpulse = delta_v_ti * delta_dirac * e_unit_dir

Fpulse(1) = delta_v_ti * delta_dirac * e_unit_dir(1)
Fpulse(2) = delta_v_ti * delta_dirac * e_unit_dir(2)
Fpulse(3) = delta_v_ti * delta_dirac * e_unit_dir(3)

!Fpulses_1 = delta_v_ti(1) * delta_dirac * er 
!Fpulses_2 = delta_v_ti(2) * delta_dirac * et
!Fpulses_3 = delta_v_ti(3) * delta_dirac * en

!Fpulses = Fpulses_1 + Fpulses_2 + Fpulses_3
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Partial derivatives w.r.t state vector
PDr = 0.0D0
PDv = 0.0D0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Partial derivatives w.r.t unknown parameters to be estimated
! ----------------------------------------------------------------------
!PD_param = delta_dirac * e_unit_dir

PD_param(1,1) = delta_dirac * e_unit_dir(1)
PD_param(2,1) = delta_dirac * e_unit_dir(2)
PD_param(3,1) = delta_dirac * e_unit_dir(3)
! ----------------------------------------------------------------------

IF (delta_dirac == 10) THEN

print *," " 
print *,"m_pd_pulses " 
print *,"delta_dirac ", delta_dirac 
print *,"delta_v_ti ", delta_v_ti 
print *,"dir ", dir 
print *,"e_unit_dir ", e_unit_dir 
print *,"Fpulse ", Fpulse 
print *,"PD_param ", PD_param
print *," " 
 

END IF

END SUBROUTINE

END Module
