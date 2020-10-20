      program main_orbit_integration019


! ----------------------------------------------------------------------
! Program:	main_orbit_integration019.f03
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
! 	Precise Orbit Determination (POD) version: Estimation of empirical forces parameters (bias, cycle-per-rev) that lead to mm-cm level orbital accuracy w.r.t. IGS precise orbits
! - Dr. Thomas Papanikolaou, 30 January 2019
! 	POD version upgrade: Ocean tides effect revision that has significant impact on longer orbit arcs e.g. 3 days 
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_orbdet
      USE m_orbext
      USE m_writearray
      USE m_statorbit
      !USE m_readmatrixfile
      USE m_writeorbit
      USE m_writesp3_hd
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: filename, EQMfname, VEQfname				
	  REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Vres, Xsigma  
      REAL (KIND = prec_d), DIMENSION(3) :: Vrms 	    
	  REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT
! ----------------------------------------------------------------------
      CHARACTER (LEN=2) :: GNSS_id
	  INTEGER (KIND = prec_int2) :: ORB_mode 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: orbconfig, orbcomp_int
      CHARACTER (LEN=100) :: EQMfname1, EQMfname2, EQMfname3
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb1, orb2, orb3
      INTEGER (KIND = prec_int8) :: sz1, sz2 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb, dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_XYZ, dorb_RTN, dorb_Kepler
      INTEGER (KIND = prec_int2) :: AllocateStatus
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: planet_matrix
      CHARACTER (LEN=3) :: sat_prn
      INTEGER (KIND = prec_int2) :: sat_vel	  

	  
! CPU Time
CALL cpu_time (CPU_t0)


! ----------------------------------------------------------------------
! Configuration files of Orbit parameterization:
!EQMfname = 'EQM.in'
!VEQfname = 'VEQ.in'
! ----------------------------------------------------------------------
!EQMfname = 'G29_EQM_2011.in'
!VEQfname = 'G29_VEQ_2011.in'
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! IGS Orbit Integration Test
! ----------------------------------------------------------------------
! G01
!EQMfname = 'IGS_test_G01_EQM_not_nst.in'
!EQMfname = 'IGS_test_G01_EQM_not.in'
!EQMfname = 'IGS_test_G01_EQM_wot.in'
!VEQfname = 'IGS_test_G01_VEQ_wot.in'

! G02
!EQMfname = 'IGS_test_G02_EQM_not_nst.in'
!EQMfname = 'IGS_test_G02_EQM_not.in'
!EQMfname = 'IGS_test_G02_EQM_wot.in'
! ----------------------------------------------------------------------
!Print *,"Orbit Parameterisation: ", EQMfname 

! VEQ configuration files
!VEQfname = 'IGS_test_G01_VEQ_not_nst.in'
!VEQfname = 'IGS_test_G01_VEQ_not.in'
!VEQfname = 'IGS_test_G01_VEQ_wot.in'

!VEQfname = 'IGS_test_G02_VEQ_not_nst.in'
!VEQfname = 'IGS_test_G02_VEQ_not.in'
!VEQfname = 'IGS_test_G02_VEQ_wot.in'
! ----------------------------------------------------------------------
 
! ----------------------------------------------------------------------
! Configuration files triples
! ----------------------------------------------------------------------
! Single or Triple
orbconfig = 1

! WOT
!EQMfname1 = 'IGS_test_G01_EQM_wot.in'
!EQMfname1 = 'IGS_test_G02_EQM_wot.in'
! NOT-NST
EQMfname1 = 'IGS_test_G01_gfm_egm2008.in'
EQMfname1 = 'IGS_test_G02_gfm_egm2008.in'


! Configuration files triples
if (orbconfig == 1) then 
! G01
EQMfname1 = 'IGS_test_G01_EQM_not_nst.in'
EQMfname2 = 'IGS_test_G01_EQM_not.in'
EQMfname3 = 'IGS_test_G01_EQM_wot.in'

else if (orbconfig == 2) then
! G02
EQMfname1 = 'IGS_test_G02_EQM_not_nst.in'
EQMfname2 = 'IGS_test_G02_EQM_not.in'
EQMfname3 = 'IGS_test_G02_EQM_wot.in'

else if (orbconfig == 3) then
! Gravity Field models comparison
EQMfname1 = 'IGS_test_G01_gfm_gocco05s.in'
EQMfname2 = 'IGS_test_G01_gfm_eigen6s2.in'
EQMfname3 = 'IGS_test_G01_gfm_egm2008.in'

else if (orbconfig == 4) then
! Gravity Field models comparison
EQMfname1 = 'IGS_test_G02_gfm_goco05s.in' 
EQMfname2 = 'IGS_test_G02_gfm_eigen6s2.in'
EQMfname3 = 'IGS_test_G02_gfm_egm2008.in'

End IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Internal Orbit comparison
! 1. orb_icrf
! 2. orb_itrf
! 3. orbext_ICRF
orbcomp_int = 1
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
EQMfname = EQMfname1
Print *,""
Print *,"Orbit Parameterisation: ", EQMfname 
! ----------------------------------------------------------------------
! Precise Orbit Determination or Orbit Prediction
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms, Xsigma)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF: "
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  Vrms
!PRINT *,"Orbit Determination: Completed"
CALL cpu_time (CPU_t1)
!PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! External Orbit Comparison (optional)
If (ORBEXT_glb > 0) Then
Call orbext(EQMfname, orb_icrf, orb_itrf, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT)
! ----------------------------------------------------------------------
PRINT *,"External Orbit comparison"
print *,"Orbit comparison: ICRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN v", stat_RTN_extC(1, 4:6)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extC(1, 4:6)

print *,"Orbit comparison: ITRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extT(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extT(1,4:6)
End If
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit for internal comparison
if (orbcomp_int == 1) then
sz1 = size(orb_icrf, DIM = 1)
sz2 = size(orb_icrf, DIM = 2)
ALLOCATE ( orb1(sz1 , sz2), STAT = AllocateStatus)
orb1 = orb_icrf

else if (orbcomp_int == 3) then
sz1 = size(orbext_ICRF, DIM = 1)
sz2 = size(orbext_ICRF, DIM = 2)
ALLOCATE ( orb1(sz1,sz2), STAT = AllocateStatus)
orb1 = orbext_ICRF

end if
! ----------------------------------------------------------------------



IF (orbconfig>0) THEN



! ----------------------------------------------------------------------
EQMfname = EQMfname2
Print *,""
Print *,"Orbit Parameterisation: ", EQMfname 
! ----------------------------------------------------------------------
! Precise Orbit Determination or Orbit Prediction
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms, Xsigma)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF: "
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  Vrms
!PRINT *,"Orbit Determination: Completed"
CALL cpu_time (CPU_t1)
!PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! External Orbit Comparison (optional)
If (ORBEXT_glb > 0) Then
Call orbext(EQMfname, orb_icrf, orb_itrf, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT)
! ----------------------------------------------------------------------
PRINT *,"External Orbit comparison"
print *,"Orbit comparison: ICRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN v", stat_RTN_extC(1, 4:6)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extC(1, 4:6)

print *,"Orbit comparison: ITRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extT(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extT(1,4:6)
End If
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit for internal comparison
if (orbcomp_int == 1) then
sz1 = size(orb_icrf, DIM = 1)
sz2 = size(orb_icrf, DIM = 2)
ALLOCATE ( orb2(sz1 , sz2), STAT = AllocateStatus)
orb2 = orb_icrf

else if (orbcomp_int == 3) then
sz1 = size(orbext_ICRF, DIM = 1)
sz2 = size(orbext_ICRF, DIM = 2)
ALLOCATE ( orb2(sz1,sz2), STAT = AllocateStatus)
orb2 = orbext_ICRF

end if
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
EQMfname = EQMfname3
Print *,""
Print *,"Orbit Parameterisation: ", EQMfname 
! ----------------------------------------------------------------------
! Precise Orbit Determination or Orbit Prediction
CALL orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms, Xsigma)
! ----------------------------------------------------------------------
print *,"Orbit residuals in ICRF: "
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  Vrms
!PRINT *,"Orbit Determination: Completed"
CALL cpu_time (CPU_t1)
!PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! External Orbit Comparison (optional)
If (ORBEXT_glb > 0) Then
Call orbext(EQMfname, orb_icrf, orb_itrf, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT)
! ----------------------------------------------------------------------
PRINT *,"External Orbit comparison"
print *,"Orbit comparison: ICRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN v", stat_RTN_extC(1, 4:6)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extC(1, 4:6)

print *,"Orbit comparison: ITRF"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extT(1, 1:3)
!WRITE (*,FMT='(A9, 3F17.9)'),"RMS Vxyz", stat_XYZ_extT(1,4:6)
End If
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit for internal comparison
if (orbcomp_int == 1) then
sz1 = size(orb_icrf, DIM = 1)
sz2 = size(orb_icrf, DIM = 2)
ALLOCATE ( orb3(sz1 , sz2), STAT = AllocateStatus)
orb3 = orb_icrf

else if (orbcomp_int == 3) then
sz1 = size(orbext_ICRF, DIM = 1)
sz2 = size(orbext_ICRF, DIM = 2)
ALLOCATE ( orb3(sz1,sz2), STAT = AllocateStatus)
orb3 = orbext_ICRF

end if
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
! Internal Orbit comparison 
! ----------------------------------------------------------------------
Print *,""
print *,"Internal Orbit comparison: ICRF"

! NOT NST - NOT
CALL statorbit (orb1, orb2, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC)
print *,"orb1-orb2 : NOT NST - NOT"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)

! NOT NST - WOT
CALL statorbit (orb1, orb3, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC)
print *,"orb1-orb3: NOT NST - WOT"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)

! NOT - WOT
CALL statorbit (orb2, orb3, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC)
print *,"orb2-orb3 : NOT - WOT"
WRITE (*,FMT='(A9, 3F17.4)'),"RMS RTN", stat_RTN_extC(1, 1:3)
WRITE (*,FMT='(A9, 3F17.4)'),"RMS XYZ",  stat_XYZ_extC(1, 1:3)
! ----------------------------------------------------------------------
filename = "OT_dorb_icrf.out"
Call writearray (dorb_icrf, filename)
filename = "OT_dorb_RTN.out"
Call writearray (dorb_RTN, filename)
filename = "OT_dorb_Kepler.out"
Call writearray (dorb_Kepler, filename)



END IF 




! ----------------------------------------------------------------------
! Write orbit matrices to ouput files (ascii)
PRINT *,"Write orbit matrices to output files"
! ----------------------------------------------------------------------
! Estimated Orbit or Predicted Orbit
filename = "orb_icrf.out"
!Call writearray (orb_icrf, filename)
Call writeorbit (orb_icrf, filename)
filename = "orb_itrf.out"
!Call writearray (orb_itrf, filename)
Call writeorbit (orb_itrf, filename)

! Variational Equations matrices
If (ESTIM_mode_glb > 0) then
filename = "VEQ_Smatrix.out"
Call writearray (veqSmatrix, filename)
filename = "VEQ_Pmatrix.out"
Call writearray (veqPmatrix, filename)
End IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Write orbit matrices to ouput files in sp3 format
PRINT *,"Write orbit matrices in sp3 format files"
! ----------------------------------------------------------------------
! Estimated Orbit or Predicted Orbit
filename = "G01.sp3"
!SUBROUTINE writesp3_hd (wrtArray, filename, sat_prn, sat_vel)
!Call writearray (orb_icrf, filename)
sat_prn = 'G01'
sat_vel = 1
Call writesp3_hd (orb1, filename,sat_prn, sat_vel)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Planets orbit files
! ----------------------------------------------------------------------
!filename = 'Venus.out'
!CALL readmatrixfile (filename, planet_matrix)
!print *,"planet_matrix(1,:)", planet_matrix(1,:)
!print *,"planet_matrix(2,:)", planet_matrix(2,:)
! ----------------------------------------------------------------------


CALL cpu_time (CPU_t1)
!PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0
WRITE (*,FMT='(A9, F9.1)'),"CPU Time (sec)", CPU_t1-CPU_t0


End Program

