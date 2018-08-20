MODULE m_orbdet


! ----------------------------------------------------------------------
! MODULE: m_orbdet.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for GNSS Orbit Determiantion
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, CRC-SI
! Created:	20 April 2018
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains
	  
	  
SUBROUTINE orbdet (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)


! ----------------------------------------------------------------------
! SUBROUTINE: m_orbdet.f03
! ----------------------------------------------------------------------
! Purpose:
!  GNSS Orbit Determination 
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
!			Geoscience Australia, CRC-SI
! Created:	20 April 2018
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_planets
      USE mdl_tides
      USE m_orbinteg
      USE m_orb_estimator
      USE m_orbC2T  
      USE m_statdelta
      USE m_statorbit
      USE m_writearray
      IMPLICIT NONE
	  
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN)  :: EQMfname, VEQfname				
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: veqSmatrix, veqPmatrix  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Vres  
      REAL (KIND = prec_d), DIMENSION(3), INTENT(OUT) :: Vrms 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: filename
      INTEGER (KIND = prec_int2) :: VEQmode 
      INTEGER (KIND = prec_int2) :: ESTmode 
      INTEGER (KIND = prec_int2) :: Niter 
      INTEGER (KIND = prec_int8) :: i 
      INTEGER (KIND = prec_int8) :: sz1, sz2 
      INTEGER (KIND = prec_int8) :: Nepochs	  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Xmatrix, Wmatrix, Amatrix  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqC, veqT  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb0, veq0, veq1  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb, dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_XYZ, dorb_RTN, dorb_Kepler
!	  REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ, stat_RTN, stat_Kepler
      REAL (KIND = prec_d), DIMENSION(:), ALLOCATABLE :: RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr 	  
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
      CHARACTER (LEN=3) :: time_sys, time 
      REAL (KIND = prec_d), DIMENSION(6) :: Zest0_icrf, Zest0_itrf, Xo_estim
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Read orbit parameterization											
! ----------------------------------------------------------------------
Call prm_main (EQMfname)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Temp																		! ----------------------------------------------------------------------
SVEC_Zo_ESTIM = SVEC_Zo
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Estimator settings :: Module mdl_param.f03 global parameters
ESTmode = ESTIM_mode_glb
Niter = ESTIM_iter_glb
! ----------------------------------------------------------------------

!PRINT *,"Data reading"
! ----------------------------------------------------------------------
! Data reading: Gravitational Effects
! ----------------------------------------------------------------------
! Earth Gravity Field model
CALL prm_gravity (EQMfname)												
! Planetary/Lunar DE data 
CALL prm_planets (EQMfname)												
! Ocean Tides model
CALL prm_ocean (EQMfname)												
! ----------------------------------------------------------------------
! Pseudo-Observations: Precise Orbit (sp3) 
CALL prm_pseudobs (EQMfname)
! ----------------------------------------------------------------------
! External Orbit comparison: Precise Orbit (sp3)
!CALL prm_orbext (EQMfname)												
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Dynamic Orbit Determination 
! ----------------------------------------------------------------------
! Orbit estimation or propagation
!ESTmode = 1

! ----------------------------------------------------------------------
If (ESTmode > 0) then
! Orbit Estimation

! Iterations number of parameter estimation algorithm
!Niter = 1

! Orbit Numerical Integration: Equation of Motion and Variational Equations
Do i = 0 , Niter

!PRINT *,"Iteration:", i

! ----------------------------------------------------------------------
! Numerical Integration: Variational Equations
! ----------------------------------------------------------------------
!PRINT *,"VEQ Integration:"
VEQmode = 1
Call orbinteg (VEQfname, VEQmode, orb0, veqSmatrix, veqPmatrix)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Numerical Integration: Equation of Motion
! ----------------------------------------------------------------------
!PRINT *,"EQM Integration:"
VEQmode = 0
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit residuals; statistics ! ICRF
!CALL statorbit (orbext_ICRF, orb_icrf, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ, stat_RTN, stat_Kepler)
Call statdelta(pseudobs_ICRF, orb_icrf, dorb_icrf, RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr)
! ----------------------------------------------------------------------
!print *,"Orbit residuals (ICRF) RMS(XYZ)", RMSdsr(1:3)


! ----------------------------------------------------------------------
! Parameter estimation: Initial Conditions and orbit parameters
! ----------------------------------------------------------------------
!Call orb_estimator(orb_icrf, veqSmatrix, veqPmatrix, orbext_ICRF, Xmatrix, Wmatrix, Amatrix)			! ----------------------------------------------------------------------
Call orb_estimator(orb_icrf, veqSmatrix, veqPmatrix, pseudobs_ICRF, Xmatrix, Wmatrix, Amatrix)			! ----------------------------------------------------------------------

filename = "Amatrix.out"
Call writearray (Amatrix, filename)
filename = "Wmatrix.out"
Call writearray (Wmatrix, filename)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Temp: to be replaced by writing prm_in files (EQM + VEQ)
! ----------------------------------------------------------------------
!print *, "Xmatrix", Xmatrix
!print *,"SVEC_Zo", SVEC_Zo
Xo_estim(1:6) = Xmatrix(1:6,1)
SVEC_Zo_ESTIM = SVEC_Zo + Xo_estim
!print *, "SVEC_Zo_ESTIM Zo+Xmatrix", SVEC_Zo_ESTIM
! ----------------------------------------------------------------------

End Do
! ----------------------------------------------------------------------

End If

! ----------------------------------------------------------------------
! Orbit Propagation
! ----------------------------------------------------------------------
VEQmode = 0
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit residuals; statistics ! ICRF
Call statdelta(pseudobs_ICRF, orb_icrf, dorb_icrf, RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr)
sz1 = size(dorb_icrf, DIM = 1)
sz2 = size(dorb_icrf, DIM = 2)
ALLOCATE (Vres(sz1,5), STAT = AllocateStatus)
Vres = dorb_icrf(1:sz1,1:5)
Vrms  = RMSdsr(1:3)
! ----------------------------------------------------------------------
!print *,"Orbit residuals opt (ICRF) RMS(XYZ)", RMSdsr(1:3)

! ----------------------------------------------------------------------
! Orbit transformation to terrestrial frame: ICRF to ITRF
! ----------------------------------------------------------------------
! Time System according to global variable TIME_Scale (Module mdl_param.f03)
time_sys = TIME_SCALE
CALL orbC2T (orb_icrf, time_sys, orb_itrf)
! ----------------------------------------------------------------------


END SUBROUTINE

End

