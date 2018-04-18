      program main_orb


! ----------------------------------------------------------------------
! Program:	main_orb.f90
! ----------------------------------------------------------------------
! Purpose:
!  Dynamic Orbit Determination of GNSS satellites 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Cooperative Research Centre for Spatial Information, Australia
! Created:	13 September 2017
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_planets
      USE mdl_tides
!      USE mdl_orbinteg
!      USE mdl_statdelta
!      USE mdl_writearray
      USE m_orbinteg
      USE m_orb_estimator
      USE m_statdelta
      USE m_statorbit
      USE m_writearray
      USE m_orbC2T  
      IMPLICIT NONE


	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: INfname, filename, EQMfname, VEQfname				
      INTEGER (KIND = prec_int2) :: VEQmode 
      INTEGER (KIND = prec_int2) :: ESTmode 
      INTEGER (KIND = prec_int2) :: Niter 
      INTEGER (KIND = prec_int8) :: i 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: sz1, sz2 
      INTEGER (KIND = prec_int8) :: Nepochs	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf, veqZ, veqP  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Xmatrix, Wmatrix, Amatrix  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqC, veqT  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb0, veq0, veq1  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb, dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_XYZ, dorb_RTN, dorb_Kepler
	  REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ, stat_RTN, stat_Kepler
      REAL (KIND = prec_d), DIMENSION(:), ALLOCATABLE :: RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr 	  
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
! ----------------------------------------------------------------------
      CHARACTER (LEN=3) :: time_sys, time 
      REAL (KIND = prec_d), DIMENSION(6) :: Zest0_icrf, Zest0_itrf, Xo_estim



Print *,"Orbit Determination"


! CPU Time
CALL cpu_time (CPU_t0)


! ----------------------------------------------------------------------
! Configuration files of Orbit parameterization:
! ----------------------------------------------------------------------
EQMfname = 'EQM.in'
VEQfname = 'VEQ.in'
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Read orbit parameterization											
! ----------------------------------------------------------------------
Call prm_main (EQMfname)

! ESTmode =
! Niter =
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Temp																		! ----------------------------------------------------------------------
SVEC_Zo_ESTIM = SVEC_Zo
print *,"main_orb.f03 | SVEC_Zo", SVEC_Zo
print *,"main_orb.f03 | SVEC_Zo_ESTIM", SVEC_Zo_ESTIM
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Data reading: Gravitational Effects
! ----------------------------------------------------------------------
PRINT *,"Data reading"
! Earth Gravity Field model
CALL prm_gravity (EQMfname)												
! Planetary/Lunar DE data 
CALL prm_planets (EQMfname)												
! Ocean Tides model
CALL prm_ocean (EQMfname)												
! ----------------------------------------------------------------------
! Pseudo-Observations: External precise orbit 
CALL prm_pseudobs (EQMfname)
! ----------------------------------------------------------------------
! External Orbit comparison
CALL prm_orbext (EQMfname)												
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Dynamic Orbit Determination 
! ----------------------------------------------------------------------
! Orbit estimation or propagation
ESTmode = 1

! ----------------------------------------------------------------------
If (ESTmode > 0) then
! Orbit Estimation

! Iterations number of parameter estimation algorithm
Niter = 1

! Orbit Numerical Integration: Equation of Motion and Variational Equations
Do i = 0 , Niter

PRINT *," "
PRINT *,"Iteration:", i

! ----------------------------------------------------------------------
! Variational Equations numerical integration
! ----------------------------------------------------------------------
VEQmode = 1
PRINT *,"VEQ Integration:"
Call orbinteg (VEQfname, VEQmode, orb0, veqZ, veqP)
PRINT *,"VEQ Integration: Completed"

! Store VEQ matrices to ouput files 
filename = "orbVEQ_icrf.out"
Call writearray (orb0, filename)
filename = "VEQ_Zmatrix.out"
Call writearray (veqZ, filename)
filename = "VEQ_Pmatrix.out"
Call writearray (veqP, filename)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Equation of Motion numerical integration
! ----------------------------------------------------------------------
VEQmode = 0
PRINT *,"EQM Integration:"
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
PRINT *,"EQM Integration: Completed"
! ----------------------------------------------------------------------


if (1<0) Then
! ----------------------------------------------------------------------
! Orbit residuals; statistics ! ICRF
!CALL statorbit (orbext_ICRF, orb_icrf, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ, stat_RTN, stat_Kepler)

!Call statdelta(pseudobs_ICRF, orb_icrf, dorb_icrf, RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr)

PRINT *," "
print *,"Orbit residuals in ICRF | Iteration:", i
print *,"RMS r     ", stat_XYZ(1, 1:3)
print *,"RMS RTN r", stat_RTN(1, 1:3)
!print *,"RMS RTN v", stat_RTN(1, 4:6)
PRINT *," "
! ----------------------------------------------------------------------
end if


! ----------------------------------------------------------------------
! Parameter estimation: Initial Conditions and orbit parameters
! ----------------------------------------------------------------------
!Call orb_estimator(orb_icrf, veqZ, veqP, orbext_ICRF, Xmatrix, Wmatrix, Amatrix)			! ----------------------------------------------------------------------
Call orb_estimator(orb_icrf, veqZ, veqP, pseudobs_ICRF, Xmatrix, Wmatrix, Amatrix)			! ----------------------------------------------------------------------

filename = "Amatrix.out"
Call writearray (Amatrix, filename)
filename = "Wmatrix.out"
Call writearray (Wmatrix, filename)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Temp: to be replaced by writing prm_in files (EQM + VEQ)
! ----------------------------------------------------------------------
print *, "Xmatrix", Xmatrix
print *,"SVEC_Zo", SVEC_Zo
Xo_estim(1:6) = Xmatrix(1:6,1)
SVEC_Zo_ESTIM = SVEC_Zo + Xo_estim
print *, "SVEC_Zo_ESTIM Zo+Xmatrix", SVEC_Zo_ESTIM
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


End Do
! ----------------------------------------------------------------------

Else

! ----------------------------------------------------------------------
! Orbit integration (Equation of Motion)
VEQmode = 0
PRINT *,"Orbit Integration :: EQM :: in-progress"
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------

End If


CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0
PRINT *,"Orbit Integration: Completed"




PRINT *,"Orbit arc (sec)", orbarc



! ----------------------------------------------------------------------
! Orbit transformation to terrestrial frame: ICRF to ITRF
! ----------------------------------------------------------------------
! Time System according to global variable TIME_Scale (Module mdl_param.f03)
time_sys = TIME_SCALE
!time_sys = 'GPS'
CALL orbC2T (orb_icrf, time_sys, orb_itrf)

! Keplerian Orbit case
!sz1 = size(orb_icrf, DIM = 1)
!sz2 = size(orb_icrf, DIM = 2)
!ALLOCATE (orb_itrf(sz1,sz2), STAT = AllocateStatus)
!orb_itrf = orb_icrf
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! External Orbit comparison
! ----------------------------------------------------------------------
PRINT *,"External Orbit: Forming orbit arrays"
! External Orbit: Read configuration file of orbit parameterization
!CALL prm_orbext (EQMfname)
PRINT *,"External Orbit: Completed"
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit comparison statistics
! ----------------------------------------------------------------------
! External orbit: orbext_ICRF, orbext_ITRF, orbext_kepler
! module mdl_statdelta.f03 : statdelta subroutine
PRINT *,"External Orbit comparison"
! ----------------------------------------------------------------------
! ICRF
CALL statorbit (orbext_ICRF, orb_icrf, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ, stat_RTN, stat_Kepler)
print *,"Orbit comparison: ICRF"
print *,"RMS RTN r ", stat_RTN(1, 1:3)
print *,"RMS RTN v ", stat_RTN(1, 4:6)
print *,"RMS r     ", stat_XYZ(1, 1:3)
print *,"RMS v     ", stat_XYZ(1, 4:6)
print *,"RMS Kepler", stat_Kepler(1, 1:3)
print *,"RMS Kepler", stat_Kepler(1, 4:6)

! ITRF
Call statdelta(orbext_ITRF, orb_itrf, dorb_itrf, RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr)
print *,"Orbit comparison: ITRF"
print *,"RMS r     ", RMSdsr(1:3)
print *,"RMS v     ", RMSdsr(4:6)
!print *,"MIN r     ", MINdsr(1:3)
!print *,"MIN v     ", MINdsr(4:6)
!print *,"MAX r     ", MAXdsr(1:3)
!print *,"MAX v     ", MAXdsr(4:6)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Write orbit matrices to ascii files
PRINT *,"Write orbits to output files"
! ----------------------------------------------------------------------
! Computed dynamic orbit (integrated)
filename = "orb_icrf.out"
Call writearray (orb_icrf, filename)
filename = "orb_itrf.out"
Call writearray (orb_itrf, filename)

! External orbit
filename = "orbext_ICRF.out"
Call writearray (orbext_ICRF, filename)
filename = "orbext_ITRF.out"
Call writearray (orbext_ITRF, filename)
! ----------------------------------------------------------------------
! Orbit comparison residuals
! ICRF
filename = "dorb_icrf.out"
Call writearray (dorb_icrf, filename)
filename = "dorb_RTN.out"
Call writearray (dorb_RTN, filename)
filename = "dorb_Kepler.out"
Call writearray (dorb_Kepler, filename)
! ITRF
filename = "dorb_itrf.out"
Call writearray (dorb_itrf, filename)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Write Cnm and Snm matrices to ascii files
! ----------------------------------------------------------------------
! Cnm
filename = "Cnm.out"
!Call writearray (GFM_Cnm, filename)  
! Snm
filename = "Snm.out"
!Call writearray (GFM_Snm, filename)  
! ----------------------------------------------------------------------






CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0



End Program

