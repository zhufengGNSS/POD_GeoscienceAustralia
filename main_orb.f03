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


Print *,"Orbit Determination"

! CPU Time
CALL cpu_time (CPU_t0)

! ----------------------------------------------------------------------
! Configuration files of Orbit parameterization:
EQMfname = 'EQM.in'
VEQfname = 'VEQ.in'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellite Orbit Determination
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF : RMS(XYZ)", Vrms
PRINT *,"Orbit Determination: Completed"
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
! Estimated Orbit 
filename = "orb_icrf.out"
Call writearray (orb_icrf, filename)
filename = "orb_itrf.out"
Call writearray (orb_itrf, filename)
! Variational Equations matrices
filename = "VEQ_Smatrix.out"
Call writearray (veqSmatrix, filename)
filename = "VEQ_Pmatrix.out"
Call writearray (veqPmatrix, filename)
! ----------------------------------------------------------------------

CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


End Program

