SUBROUTINE eclipse_integstep (EQMfname, VEQfname, mjd, r_sat, v_sat, integstep_flag, integstep_initial, integstep_reduced)


! ----------------------------------------------------------------------
! SUBROUTINE: eclipse_integstep.f03
! ----------------------------------------------------------------------
! Purpose:
!  Reduce the step of the numerical integration method applied for orbit 
!  determination during eclipse seasons 
! ----------------------------------------------------------------------
! Input arguments:
! - EQMfname: 	Input configuration file name for the orbit parameterization 
! - VEQfname: 	Input configuration file name for the orbit parameterization 
! - mjd:		Modified Julian Day (MJD) in Terrestrial Time (including the fraction of the day)
! - r_sat: 		Satellite position vector (m) in ICRF
! - v_sat: 		Satellite velocity vector (m/sec) in ICRF
!
! Output arguments:
! - integstep_flag: 	Flag regarding the change of the orbit integration step 
! - integstep_initial:	Initial value of the orbit integration method 
! - integstep_reduced:	Reduced value of the orbit integration method 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
! 			Geoscience Australia, Frontier-SI
! Created:	11 June 2019
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN)  :: EQMfname, VEQfname				
      REAL (KIND = prec_d), INTENT(IN) :: mjd, r_sat(3), v_sat(3)!, r_sun(3)
! ----------------------------------------------------------------------
! OUT
      LOGICAL, INTENT(OUT) :: integstep_flag
      REAL (KIND = prec_d), INTENT(OUT) :: integstep_initial, integstep_reduced
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: rbody(3)
      REAL (KIND = prec_d) :: r_sun(3)
      REAL (KIND = prec_d) :: beta
      DOUBLE PRECISION  JD, Zbody(6)
      INTEGER  NTARG, NCTR, NTARG_body
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
! ----------------------------------------------------------------------


integstep_flag = .FALSE.

! ----------------------------------------------------------------------
! Julian Day Number of the input epoch
JD = mjd + 2400000.5D0
! Center celestial body: Earth
NCTR = 3 
! Sun (NTARG) Cartesian coordinates w.r.t. Earth (NCTR)
NTARG = 11
CALL  PLEPH ( JD, NTARG, NCTR, Zbody )
! Cartesian coordinates of the celestial body in meters: KM to M
rbody(1) = Zbody(1) * 1000.D0
rbody(2) = Zbody(2) * 1000.D0
rbody(3) = Zbody(3) * 1000.D0
r_sun = rbody
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Beta angle (degrees)
CALL beta_angle (r_sat, v_sat, r_sun, beta)
! ----------------------------------------------------------------------
print *,"beta", beta

! ----------------------------------------------------------------------
! Criteria for changing orbit integration stepsize
IF ( abs(beta) <= 14) THEN
	integstep_initial = integstep
IF (integstep_initial >= 300.D0) THEN
	integstep_reduced = 100.0D0
	integstep_flag = .TRUE.
END IF

write (fname_id, *) '_INT'
param_id = 'integrator_step'
write (param_value, *) integstep_reduced
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

END IF 
! ----------------------------------------------------------------------

END SUBROUTINE
