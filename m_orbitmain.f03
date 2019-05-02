MODULE m_orbitmain


! ----------------------------------------------------------------------
! MODULE: m_orbitmain.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for Precise Orbit Determination
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	21 March 2019
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains
	  
	  
SUBROUTINE orbitmain (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)

! ----------------------------------------------------------------------
! SUBROUTINE:	orbitmain.f03
! ----------------------------------------------------------------------
! Purpose:
!  Precise Orbit Determination 
! ----------------------------------------------------------------------
! Input arguments:
! - EQMfname: 	Input cofiguration file name for the orbit parameterization 
! - VEQfname: 	Input cofiguration file name for the orbit parameterization 
!
! Output arguments:
! - orb_icrf: 	Satellite orbit array in ICRF including the following per epoch:
!               - Modified Julian Day number (including the fraction of the day) 
!				- Seconds since 00h 
!				- Position vector (m)
!				- Velocity vector (m/sec)
! - orb_itrf: 	Satellite orbit array in ITRF including the following per epoch:
!               - Modified Julian Day number (including the fraction of the day) 
!				- Seconds since 00h 
!				- Position vector (m)
!				- Velocity vector (m/sec)
! - veqSmatrix:	State trasnition matrix obtained from the Variational Equations solution based on numerical integration methods
! - veqPmatrix: Sensitivity matrix obtained from the Variational Equations solution based on numerical integration methods
! ----------------------------------------------------------------------
! Note 1:
! The time scale of the 2 first collumns of the orbit arrays (MJD and Seoncds since 00h) 
! refer to the time system defined by the global variable TIME_SCALE in the module mdl_param.f03
! according to the input parameterization file 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	21 March 2019
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_orbdet
      USE m_orbext
      USE m_writearray
!      USE m_writeorbit
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: filename, EQMfname, VEQfname				
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Vres  
      REAL (KIND = prec_d), DIMENSION(3) :: Vrms 	    
	  REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT
! ----------------------------------------------------------------------
      CHARACTER (LEN=2) :: GNSS_id
	  INTEGER (KIND = prec_int2) :: ORB_mode
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Nsat, isat
	  INTEGER (KIND = prec_int8) :: iepoch, iparam
	  INTEGER (KIND = prec_int8) :: i
	  INTEGER (KIND = prec_int8) :: sz1, sz2, Nepochs, N2_orb, N2_veqSmatrix, N2_veqPmatrix, N2sum  
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbit_veq  
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus  
	  CHARACTER (LEN=3), ALLOCATABLE :: PRN_array(:)
	  CHARACTER (LEN=3) :: PRN_isat
	  INTEGER :: ios
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Precise Orbit Determination or Orbit Prediction
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF:" 
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ", Vrms
!PRINT *,"Orbit Determination: Completed"
!CALL cpu_time (CPU_t1)
!PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


If (ORBEXT_glb > 0) Then
! ----------------------------------------------------------------------
! External Orbit Comparison (optional)
Call orbext(EQMfname, orb_icrf, orb_itrf, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT)
! ----------------------------------------------------------------------
PRINT *,"External Orbit comparison"
print *,"Orbit comparison: ICRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extC(1, 4:6)

print *,"Orbit comparison: ITRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extT(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extT(1,4:6)
End If


! ----------------------------------------------------------------------
! Write orbit matrices to output files (ascii)
!PRINT *,"Write orbit matrices to output files"
! ----------------------------------------------------------------------
! Estimated Orbit or Predicted Orbit
filename = "orb_icrf.out"
!Call writearray (orb_icrf, filename)
!Call writeorbit (orb_icrf, filename)
filename = "orb_itrf.out"
!Call writearray (orb_itrf, filename)
!Call writeorbit (orb_itrf, filename)

! Variational Equations matrices
If (ESTIM_mode_glb > 0) then
filename = "VEQ_Smatrix.out"
Call writearray (veqSmatrix, filename)
filename = "VEQ_Pmatrix.out"
Call writearray (veqPmatrix, filename)
End IF
! ----------------------------------------------------------------------



End SUBROUTINE


End MODULE
