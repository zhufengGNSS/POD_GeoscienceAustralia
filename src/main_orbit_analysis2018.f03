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
! Last modified: 
! - Dr. Thomas Papanikolaou, 3 May 2018
! 	Preliminary version of GNSS dynamic orbit determination	
! - Dr. Thomas Papanikolaou, 25 June 2018
! 	Version with minor revisions
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

! ----------------------------------------------------------------------
! Orbit analysis IGS2018
! Configuration files:
! ----------------------------------------------------------------------
! GNSS constellations option
! ----------------------------------------------------------------------
! IGS Final Precise Orbits
GNSS_id = 'GI'
! MGEX
!GNSS_id = 'G'
!GNSS_id = 'R'
!GNSS_id = 'E'
!GNSS_id = 'CI' ! IGSO
!GNSS_id = 'CM' ! MEO
!GNSS_id = 'CG' ! GEO
! ----------------------------------------------------------------------
! Orbit Determination or Prediction (1/0)
ORB_mode = 1
! ----------------------------------------------------------------------
If (GNSS_id == 'GI') Then
! GPS
if (ORB_mode == 1) Then
! OD
EQMfname = 'G29_EQM_2011.in'
VEQfname = 'G29_VEQ_2011.in'
!EQMfname = 'G29_EQM_igs.in'
!VEQfname = 'G29_VEQ_igs.in'
Else
! OP
EQMfname = 'G29_OP_2011.in'
VEQfname = 'VEQ.in'
end if


Else If (GNSS_id == 'G') Then
! GPS
if (ORB_mode == 1) Then
! OD
EQMfname = 'G29_EQM.in'
VEQfname = 'G29_VEQ.in'
Else
! OP
EQMfname = 'G29_OP.in'
VEQfname = 'VEQ.in'
end if


Else If (GNSS_id == 'R') Then
! GLONASS
if (ORB_mode == 1) Then
! OD
EQMfname = 'R07_EQM.in'
VEQfname = 'R07_VEQ.in'
Else
! OP
EQMfname = 'R07_OP.in'
VEQfname = 'VEQ.in'
End if


Else If (GNSS_id == 'E') Then
! Galileo
if (ORB_mode == 1) Then
! OD
EQMfname = 'E19_EQM.in'
VEQfname = 'E19_VEQ.in'
Else
! OP
EQMfname = 'E19_OP.in'
VEQfname = 'VEQ.in'
end if


Else If (GNSS_id == 'CI') Then
! Beidou IGSO
if (ORB_mode == 1) Then
! OD
!EQMfname = 'IGSO_C07_EQM.in'
!VEQfname = 'IGSO_C07_VEQ.in'
EQMfname = 'C07_EQM_IGSO.in'
VEQfname = 'C07_VEQ_IGSO.in'
Else
! OP
EQMfname = 'C07_OP_IGSO.in'
VEQfname = 'VEQ.in'
End if


Else If (GNSS_id == 'CM') Then
! Beidou MEO
if (ORB_mode == 1) Then
! OD
EQMfname = 'C11_EQM_MEO.in'
VEQfname = 'C11_VEQ_MEO.in'
Else
! OP
EQMfname = 'C11_OP_MEO.in'
VEQfname = 'VEQ.in'
End if


Else If (GNSS_id == 'CG') Then
! Beidou GEO
if (ORB_mode == 1) Then
! OD
EQMfname = 'C03_EQM_GEO.in'
VEQfname = 'C03_VEQ_GEO.in'
Else
! OP
EQMfname = 'C03_OP_GEO.in'
VEQfname = 'VEQ.in'
End if


End if

print *,"Orbit configuration files ", EQMfname
print *,"Orbit configuration files ", VEQfname
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Precise Orbit Determination or Orbit Prediction
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF : RMS(XYZ)", Vrms
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
print *,"RMS RTN", stat_RTN_extC(1, 1:3)
!print *,"RMS RTN v", stat_RTN_extC(1, 4:6)
print *,"RMS XYZ", stat_XYZ_extC(1, 1:3)
!print *,"RMS Vxyz", stat_XYZ_extC(1, 4:6)

print *,"Orbit comparison: ITRF"
print *,"RMS XYZ", stat_XYZ_extT(1,1:3)
!print *,"RMS Vxyz", stat_XYZ_extT(1,4:6)
End If


! ----------------------------------------------------------------------
! Write orbit matrices to ouput files (ascii)
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

