      program main_orb


! ----------------------------------------------------------------------
! Program:	main_orb.f90
! ----------------------------------------------------------------------
! Purpose:
!  Dynamic Orbit Determination of GNSS satellites 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, CRC-SI
! Created:	13 September 2017
! ----------------------------------------------------------------------
! POD version major modifications highlights: 
! Last modified 
! - Dr. Thomas Papanikolaou, 3 May 2018
! 	Preliminary version of GNSS dynamic orbit determination	
! - Dr. Thomas Papanikolaou, 25 June 2018
! 	Version with minor revisions
! - Dr. Thomas Papanikolaou, 30 November 2018
! 	Precise Orbit Determination (POD) version including the estimation of empirical parameters (forces related)
! - Dr. Thomas Papanikolaou, 30 January 2019
! 	POD version including the upgrade of the ocean tides effect with impact on longer orbit arcs e.g. 3 days 
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_orbdet
      USE m_orbext
      USE m_writearray
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

	  
! CPU Time
CALL cpu_time (CPU_t0)


! ----------------------------------------------------------------------
! Configuration files of Orbit parameterization:
EQMfname = 'EQM.in'
VEQfname = 'VEQ.in'
! ----------------------------------------------------------------------
!EQMfname = 'G29_EQM_2011.in'
!VEQfname = 'G29_VEQ_2011.in'
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Precise Orbit Determination or Orbit Prediction
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF:" 
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ", Vrms
!PRINT *,"Orbit Determination: Completed"
CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


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
PRINT *,"Write orbit matrices to output files"
! ----------------------------------------------------------------------
! Estimated Orbit or Predicted Orbit
filename = "orb_icrf.out"
Call writearray (orb_icrf, filename)
filename = "orb_itrf.out"
Call writearray (orb_itrf, filename)

! Variational Equations matrices
If (ESTIM_mode_glb > 0) then
filename = "VEQ_Smatrix.out"
Call writearray (veqSmatrix, filename)
filename = "VEQ_Pmatrix.out"
Call writearray (veqPmatrix, filename)
End IF
! ----------------------------------------------------------------------




CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


End Program

