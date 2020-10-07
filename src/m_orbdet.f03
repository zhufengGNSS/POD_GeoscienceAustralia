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
	  
	  
SUBROUTINE orbdet (EQMfname, VEQfname, orb_icrf_final, orb_itrf_final, veqSmatrix_final, veqPmatrix_final, Vres, Vrms)


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
! - Vres:		Orbit residuals matrix
! - Vrms: 		RMS values of the orbit residuals per X,Y,Z components 
! ----------------------------------------------------------------------
! Note 1:
! The time scale of the 2 first collumns of the orbit arrays (MJD and Seoncds since 00h) 
! refer to the time system defined by the global variable TIME_SCALE in the module mdl_param.f03
! according to the input parameterization file 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, CRC-SI
! Created:	20 April 2018
!
! Changes:  18-12-2018  Tzupang Tseng : Enabled the function of the ECOM SRP estimation and added some conditions to judge which model
!                                       is used to improve the GNSS orbit modelling (Currently the ECOM model is only estimated with full
!                                       9 coefficients or 3 bias terms. The adjustable function has not been ready yet)
!           21-02-2019  Tzupang Tseng : The adjustable function of the ECOM model has been activated.
!           06-08-2019  Tzupang Tseng : Added a function to skip bad orbits with zero value in SP3 file
!           03-12-2019  Tzupang Tseng : Added a function of estimating parameters in simple box wing model
!
! Last modified:
! 20 May 2019,	Dr. Thomas Papanikolaou
!				General upgrade for supporting the POD Tool modes of orbit determination and prediction 
! 				1. Orbit Determination (pseudo-observations; orbit fitting)
! 				2. Orbit Determination and Prediction
! 				3. Orbit Integration (Equation of Motion only)
! 				4. Orbit Integration and Partials (Equation of Motion and Variational Equations)
! 30 May 2019,	Dr. Thomas Papanikolaou
! 				Minor modification for the case of Orbit Determination and Prediction (POD mode 2) 
!				for computing the orbit residuals of the estimated part only of the orbits without the predicted part 
! 11 June 2019,	Dr. Thomas Papanikolaou
! 				Modification of the orbit integrator step in case of eclipse seasons; 
!               Resize the orbit and partials matrices according to the initial integrator step prior passing to output arguments 
! 17/08/2020,	Dr. Thomas Papanikolaou, Pseudo-stochastic pulses (velocity) added  
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_config
      USE mdl_planets
      USE mdl_tides
      USE m_orbinteg
      USE m_orb_estimator
      USE m_orbC2T  
      USE m_statdelta
      USE m_statorbit
      USE m_writearray
      USE m_orbresize
      USE m_matrixreverse
      USE m_matrixmerge
      USE m_ecom_init
      IMPLICIT NONE
	  
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN)  :: EQMfname, VEQfname				
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orb_icrf_final, orb_itrf_final  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: veqSmatrix_final, veqPmatrix_final  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Vres  
      REAL (KIND = prec_d), DIMENSION(3), INTENT(OUT) :: Vrms 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------  
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: filename
      CHARACTER (LEN=10) :: DOYSTR
      INTEGER (KIND = prec_int2) :: VEQmode 
      INTEGER (KIND = prec_int2) :: ESTmode 
      INTEGER (KIND = prec_int2) :: Niter,srp_i 
      INTEGER (KIND = prec_int8) :: i, j, k, ii, PD_Param_ID 
      INTEGER (KIND = prec_int8) :: sz1, sz2 
      INTEGER (KIND = prec_int8) :: Nepochs	  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Xmatrix, Wmatrix, Amatrix  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqC, veqT  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb0, veq0, veq1  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb, dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_XYZ, dorb_RTN, dorb_Kepler
	  REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ, stat_RTN, stat_Kepler
      REAL (KIND = prec_d), DIMENSION(:), ALLOCATABLE :: RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr 	  
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
      CHARACTER (LEN=3) :: time_sys, time 
      REAL (KIND = prec_d), DIMENSION(6) :: Zest0_icrf, Zest0_itrf, Xo_estim
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: Bias_corr(3), CPR_corr(3,2)
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: fname, fname1, fname2				
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
      REAL (KIND = prec_d) :: apriori_3(3), apriori_2(2) 				
	  REAL (KIND = prec_q) :: Bias_0(3)
	  REAL (KIND = prec_q) :: CPR_CS_0(3,2)
      REAL (KIND = prec_q), DIMENSION(:), ALLOCATABLE :: ECOM_0_coef,ECOM_coef
      REAL (KIND = prec_q), DIMENSION(:), ALLOCATABLE :: ECOM_accel_aposteriori
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: EQMfname_pred, VEQfname_pred
      REAL (KIND = prec_d) :: orbarc_sum
      INTEGER (KIND = prec_int8) :: Nepochs_estim	  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf_estim
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: mjd, r_sat(3), v_sat(3)
      LOGICAL :: integstep_flag
      REAL (KIND = prec_d) :: integstep_initial, integstep_reduced, integstep_orb, integstep_orbback
      INTEGER (KIND = prec_int8) :: integstep_rate 
      INTEGER (KIND = prec_int8) :: Nepochs_0, Nepochs_stepsmall, n2_orb, n2_veqs, n2_veqp 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_stepsmall, veqSmatrix_stepsmall, veqPmatrix_stepsmall
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix  
      INTEGER (KIND = prec_int2) :: VEQ_refsys
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: orbintegr_back_flag
      CHARACTER (LEN=100) :: EQMfname_back, VEQfname_back
      REAL (KIND = prec_d) :: orbarc_back
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_back, veqSmatrix_back, veqPmatrix_back
      CHARACTER (LEN=100):: mesg
      INTEGER (KIND = prec_int2) :: pseudobs_opt
! ----------------------------------------------------------------------
	  REAL (KIND = prec_q) :: delta_v_corr(3)
      INTEGER (KIND = prec_int8) :: PULSE_param, N_PULSE_param
      INTEGER (KIND = prec_int8) :: i_pulse
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Variable initialisation
CPR_corr = 0.d0
Bias_corr = 0.d0

! ----------------------------------------------------------------------
! TEMP: Pseudo-stochastic pulses
N_PULSE_param_glb = 3
PULSE_param_glb = 1
MJD_DELTA_V_glb = MJD_to + 0.5D0
!DELTA_V_apriori_glb = (/ 0.0D0, 0.0D0, 0.0D0 /) 
DELTA_V_aposteriori_glb = (/ 1.0D-9, 1.0D-9, 1.0D-9 /) 
!DELTA_V_aposteriori_glb = (/ 0.0D0, 0.0D0, 0.0D0 /) 
! ----------------------------------------------------------------------
PULSE_param = PULSE_param_glb
	  
! ----------------------------------------------------------------------
! Read orbit parameterization											
! ----------------------------------------------------------------------
Call prm_main (EQMfname)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Temp																		! ----------------------------------------------------------------------
SVEC_Zo_ESTIM = SVEC_Zo
!Bias_accel_aposteriori = Bias_accel_glb
!CPR_CS_aposteriori = CPR_CS_glb
! ----------------------------------------------------------------------	

! ----------------------------------------------------------------------
! Estimator settings :: Module mdl_param.f03 global parameters
ESTmode = ESTIM_mode_glb
Niter = ESTIM_iter_glb
! ----------------------------------------------------------------------

!Print *,"Orbit ESTmode:", ESTmode
If (ESTmode == 0) then
Print *,"Orbit Propagation"
Else 
Print *,"Orbit Determination"
End IF

!PRINT *,"Data reading"
! ----------------------------------------------------------------------
! Data reading: Gravitational Effects
! ----------------------------------------------------------------------
! Earth Gravity Field model
!CALL prm_gravity (EQMfname)												
! Planetary/Lunar DE data 
!CALL prm_planets (EQMfname)												
! Ocean Tides model
!CALL prm_ocean (EQMfname)												
! ----------------------------------------------------------------------
! Pseudo-Observations: Precise Orbit (sp3) 
pseudobs_opt = 1
CALL prm_pseudobs (EQMfname, pseudobs_opt)
! ----------------------------------------------------------------------
! External Orbit comparison: Precise Orbit (sp3)
!CALL prm_orbext (EQMfname)												
! ----------------------------------------------------------------------
! Skip bad orbits with zero value in SP3 file
!CALL scan0orb
! ----------------------------------------------------------------------
! Reference system of Variational Equations solution' matrices (Smatrix, Pmatrix)
! and orbit parameter estimation 
! ----------------------------------------------------------------------
VEQ_refsys = 0
IF (VEQ_REFSYS_cfg == 'ICRS') THEN 
	VEQ_refsys = 1 
ELSE IF (VEQ_REFSYS_cfg == 'ITRS') THEN
	VEQ_refsys = 2 
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!  Control of orbit integrator step during eclipse seasons 
! ----------------------------------------------------------------------
! Global variables via mdl_param.f90
mjd = MJD_to
r_sat(1:3) = SVEC_Zo(1:3)
v_sat(1:3) = SVEC_Zo(4:6)
CALL eclipse_integstep (EQMfname, VEQfname, mjd, r_sat, v_sat, integstep_flag, integstep_initial, integstep_reduced)
! ----------------------------------------------------------------------
!print *,"integstep_flag,integstep_initial, integstep_reduced", integstep_flag,integstep_initial, integstep_reduced 

! ----------------------------------------------------------------------
! Initial conditions
! ----------------------------------------------------------------------
! Empirical parameters apriori values set to zero
IF (EMP_param_glb > 0) THEN
   i = 999
   Bias_0 = (/ 0.0D0, 0.0D0, 0.0D0/)
   CPR_CS_0(1,:) = (/ 0.0D0, 0.0D0/)
   CPR_CS_0(2,:) = (/ 0.0D0, 0.0D0/)
   CPR_CS_0(3,:) = (/ 0.0D0, 0.0D0/)
!Call empirical_init (i, Bias_0, CPR_CS_0)
   Call empirical_init_file (i, Bias_0, CPR_CS_0)
END IF
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Initial conditions for solar radiation pressure
! ----------------------------------------------------------------------
IF(ECOM_param_glb == 1) PRINT*,'ECOM1 SRP MODEL IS ACTIVATED'
IF(ECOM_param_glb == 2) PRINT*,'ECOM2 SRP MODEL IS ACTIVATED'
IF(ECOM_param_glb == 3) PRINT*,'SIMPLE BOX WING IS ACTIVATED'
IF(ECOM_param_glb > 3)  PRINT*,'UNKNOWN SRP MODEL IS ACTIVATED :-('
!   ALLOCATE (ECOM_coef(NPARAM_glb), STAT = AllocateStatus)
   ALLOCATE (ECOM_coef(ECOMNUM), STAT = AllocateStatus)
   if (AllocateStatus .ne. 0) then
        write(mesg, *) "Not enough memory - failed to allocate ECOM_coef, dimension = ", NPARAM_EMP_ECOM_glb
        call report('FATAL', pgrm_name, 'orbdet', mesg, 'src/m_orbdet.f03', 1)
   end if
!   ALLOCATE (ECOM_accel_aposteriori(NPARAM_glb), STAT = AllocateStatus)
   ALLOCATE (ECOM_accel_aposteriori(ECOMNUM), STAT = AllocateStatus)
   if (AllocateStatus .ne. 0) then
        write(mesg, *) "Not enough memory - failed to allocate ECOM_accel_aposteroiri, ", &
                "dimension = ", NPARAM_EMP_ECOM_glb
        call report('FATAL', pgrm_name, 'orbdet', mesg, 'src/m_orbdet.f03', 1)
   end if
   CALL ecom_init (ECOM_0_coef)

IF(EMP_param_glb /= 0) PRINT*,'EMPIRICAL MODEL IS ACTIVATED '
IF(EMP_param_glb==0 .AND. ECOM_param_glb==0)PRINT*,'NEITHER ECOM SRP MODEL or EMPIRICAL MODEL ARE ACTIVATED'

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Dynamic Orbit Determination 
! ----------------------------------------------------------------------
! POD modes :: 1, 2
! 1. Orbit Determination (pseudo-observations; orbit fitting)
! 2. Orbit Determination and Prediction
! ----------------------------------------------------------------------
If (ESTmode > 0) then
! Orbit Estimation

! Iterations number of parameter estimation algorithm
Do i = 0 , Niter
!PRINT *,"Iteration:", i

! ----------------------------------------------------------------------
! A-priori values of Parameters
! ----------------------------------------------------------------------
! Pseudo-stachastic pulses
DELTA_V_apriori_glb = DELTA_V_aposteriori_glb
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Numerical Integration: Equation of Motion and Variational Equations
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
IF (VEQ_refsys == 1) THEN
	! Orbit parameter estimator
	Call orb_estimator(orb_icrf, veqSmatrix, veqPmatrix, pseudobs_ICRF, Xmatrix, Wmatrix, Amatrix)			! ----------------------------------------------------------------------
ELSE IF (VEQ_refsys == 2) THEN
	! Time System according to global variable TIME_Scale (Module mdl_param.f03)
	time_sys = TIME_SCALE
	! Orbit transformation to terrestrial frame: ICRF to ITRF
	CALL orbC2T (orb_icrf, time_sys, orb_itrf)
	! Orbit parameter estimator
	Call orb_estimator(orb_itrf, veqSmatrix, veqPmatrix, pseudobs_ITRF, Xmatrix, Wmatrix, Amatrix)		
END IF

!filename = "Amatrix.out"
!Call writearray (Amatrix, filename)
!filename = "Wmatrix.out"
!Call writearray (Wmatrix, filename)
!filename = "pseudobs_ITRF.out"
!Call writearray (pseudobs_ITRF, filename)
!filename = "pseudobs_ICRF.out"
!Call writearray (pseudobs_ICRF, filename)
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! Temp: to be replaced by writing prm_in files (EQM + VEQ)								! ----------------------------------------------------------------------
! ----------------------------------------------------------------------   
!print *, "Xmatrix Z", Xmatrix(1:6,1)
!print *, "Xmatrix P", Xmatrix(7:NPARAM_glb+6,1)
!print *,"SVEC_Zo", SVEC_Zo
Xo_estim(1:6) = Xmatrix(1:6,1)
SVEC_Zo_ESTIM = SVEC_Zo + Xo_estim
!print *, "SVEC_Zo_ESTIM Zo+Xmatrix", SVEC_Zo_ESTIM


! ----------------------------------------------------------------------
! Empirical model
! **********************************************************************
!If (NPARAM_glb /=0) Then !(remarked by Dr. Tzupang Tseng 11-12-2018)
If (EMP_param_glb == 1 .and. ECOM_param_glb == 0) Then
        If (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 1) Then
        ! Bias & CPR
        Bias_corr = Xmatrix(7:9,1)
        CPR_corr(1,1) = Xmatrix(10,1)
        CPR_corr(1,2) = Xmatrix(11,1)
        CPR_corr(2,1) = Xmatrix(12,1)
        CPR_corr(2,2) = Xmatrix(13,1)
        CPR_corr(3,1) = Xmatrix(14,1)
        CPR_corr(3,2) = Xmatrix(15,1)
        Else If (EMP_Bias_glb(1) == 0 .and. EMP_CPR_glb(1) == 1) Then
        ! CPR
        CPR_corr(1,1) = Xmatrix(7,1)
        CPR_corr(1,2) = Xmatrix(8,1)
        CPR_corr(2,1) = Xmatrix(9,1)
        CPR_corr(2,2) = Xmatrix(10,1)
        CPR_corr(3,1) = Xmatrix(11,1)
        CPR_corr(3,2) = Xmatrix(12,1)
        Else If (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 0) Then
        ! Bias 
        Bias_corr = Xmatrix(7:9,1)
        End If

Bias_accel_aposteriori = Bias_accel_glb + Bias_corr

!print *, "Bias_accel_aposteriori", Bias_accel_aposteriori
!print *, "Bias_accel_glb", Bias_accel_glb
!print *, "Bias_corr", Bias_corr

CPR_CS_aposteriori = CPR_CS_glb + CPR_corr
!print *, "CPR_CS_aposteriori", CPR_CS_aposteriori
!print *, "CPR_CS_glb", CPR_CS_glb
!print *, "CPR_corr", CPR_corr

! ----------------------------------------------------------------------
! Empirical parameters
! Bias and CPR terms
Bias_0 = Bias_accel_aposteriori
CPR_CS_0 = CPR_CS_aposteriori
Call empirical_init (i, Bias_0, CPR_CS_0)
End If  ! End of empirical model
! **********************************************************************
! ----------------------------------------------------------------------
! ECOM-based SRP model
! **********************************************************************
! added by Dr. Tzupang Tseng 11-12-2018
If (ECOM_param_glb == 1 .or. ECOM_param_glb == 2 .and. EMP_param_glb == 0) Then
        PD_Param_ID = 0
        If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        End IF
        If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        End IF
        If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        End IF
        If (ECOM_CPR_glb(1) == 1) THEN
        ! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        ! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        End IF
        If (ECOM_CPR_glb(2) == 1) THEN
        ! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        ! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        End IF
        If (ECOM_CPR_glb(3) == 1) THEN
        ! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        ! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_coef (PD_Param_ID) = Xmatrix(6+PD_Param_ID,1)
        End If

        IF (ECOMNUM /= PD_Param_ID) THEN
        PRINT*, 'THE NUMBER OF ECOM PARAMETERS IS NOT CONSISTENT'
        PRINT*,           'ECOMNUM  =', ECOMNUM
        PRINT*,           'PD_Param_ID =', PD_Param_ID
        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
        STOP
        END IF
END IF

IF (ECOM_param_glb == 3 .and. EMP_param_glb == 0) THEN
        PD_Param_ID = 7   
        DO k = 1,PD_Param_ID
        ECOM_coef (k) = Xmatrix(6+k,1)
        END DO

        IF (ECOMNUM /= PD_Param_ID) THEN
        PRINT*, 'THE NUMBER OF ECOM PARAMETERS IS NOT CONSISTENT'
        PRINT*,           'ECOMNUM  =', ECOMNUM
        PRINT*,           'PD_Param_ID =', PD_Param_ID
        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
        STOP
        END IF
END IF

! Switch on both ECOM and empirical models
! -------------------------------------------------------------------
IF (ECOM_param_glb /= 0 .AND. EMP_param_glb == 1) THEN
        If (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 1) Then
        ! Bias & CPR
        Bias_corr = Xmatrix(7:9,1)
        CPR_corr(1,1) = Xmatrix(10,1)
        CPR_corr(1,2) = Xmatrix(11,1)
        CPR_corr(2,1) = Xmatrix(12,1)
        CPR_corr(2,2) = Xmatrix(13,1)
        CPR_corr(3,1) = Xmatrix(14,1)
        CPR_corr(3,2) = Xmatrix(15,1)
        Else If (EMP_Bias_glb(1) == 0 .and. EMP_CPR_glb(1) == 1) Then
        ! CPR
        CPR_corr(1,1) = Xmatrix(7,1)
        CPR_corr(1,2) = Xmatrix(8,1)
        CPR_corr(2,1) = Xmatrix(9,1)
        CPR_corr(2,2) = Xmatrix(10,1)
        CPR_corr(3,1) = Xmatrix(11,1)
        CPR_corr(3,2) = Xmatrix(12,1)
        Else If (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 0) Then
        ! Bias
        Bias_corr = Xmatrix(7:9,1)
        End If
        Bias_accel_aposteriori = Bias_accel_glb + Bias_corr
        CPR_CS_aposteriori = CPR_CS_glb + CPR_corr
        Bias_0 = Bias_accel_aposteriori
        CPR_CS_0 = CPR_CS_aposteriori
        Call empirical_init (i, Bias_0, CPR_CS_0)

        If (ECOM_param_glb == 1 .or. ECOM_param_glb == 2) Then
        PD_Param_ID = 0
                IF (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 1) THEN
                        If (ECOM_Bias_glb(1) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        End IF
                        If (ECOM_Bias_glb(2) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        End IF
                        If (ECOM_Bias_glb(3) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(1) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(2) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(3) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(15+PD_Param_ID,1)
                        End If

                        IF (NPARAM_glb /= PD_Param_ID + 9) THEN
                        PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
                        PRINT*,           'NPARAM_glb  =', NPARAM_glb
                        PRINT*,           'PD_Param_ID =', PD_Param_ID+9
                        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
                        STOP
                        END IF

                ELSE IF (EMP_Bias_glb(1) == 0 .and. EMP_CPR_glb(1) == 1) THEN
                        If (ECOM_Bias_glb(1) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        End IF
                        If (ECOM_Bias_glb(2) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        End IF
                        If (ECOM_Bias_glb(3) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(1) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(2) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(3) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(12+PD_Param_ID,1)
                        End If

                        IF (NPARAM_glb /= PD_Param_ID + 6) THEN
                        PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
                        PRINT*,           'NPARAM_glb  =', NPARAM_glb
                        PRINT*,           'PD_Param_ID =', PD_Param_ID+6
                        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
                        STOP
                        END IF

                ELSE IF (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 0) THEN
                        If (ECOM_Bias_glb(1) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        End IF
                        If (ECOM_Bias_glb(2) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        End IF
                        If (ECOM_Bias_glb(3) == 1) Then
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(1) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(2) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        End IF
                        If (ECOM_CPR_glb(3) == 1) THEN
                        ! C term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        ! S term
                        PD_Param_ID = PD_Param_ID + 1
                        ECOM_coef (PD_Param_ID) = Xmatrix(9+PD_Param_ID,1)
                        End If

                        IF (NPARAM_glb /= PD_Param_ID + 3) THEN
                        PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
                        PRINT*,           'NPARAM_glb  =', NPARAM_glb
                        PRINT*,           'PD_Param_ID =', PD_Param_ID+3
                        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
                        STOP
                        END IF

                END IF
        END IF

        IF (ECOM_param_glb == 3) THEN
        PD_Param_ID = 7
                IF (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 1) THEN
                        DO k = 1,PD_Param_ID
                        ECOM_coef (k) = Xmatrix(15+k,1)
                        END DO

                        IF (NPARAM_glb /= PD_Param_ID + 9) THEN
                        PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
                        PRINT*,           'NPARAM_glb  =', NPARAM_glb
                        PRINT*,           'PD_Param_ID =', PD_Param_ID+9
                        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
                        STOP
                        END IF

                ELSE IF (EMP_Bias_glb(1) == 0 .and. EMP_CPR_glb(1) == 1) THEN
                        DO k = 1,PD_Param_ID
                        ECOM_coef (k) = Xmatrix(12+k,1)
                        END DO

                        IF (NPARAM_glb /= PD_Param_ID + 6) THEN
                        PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
                        PRINT*,           'NPARAM_glb  =', NPARAM_glb
                        PRINT*,           'PD_Param_ID =', PD_Param_ID+6
                        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
                        STOP
                        END IF

                ELSE IF (EMP_Bias_glb(1) == 1 .and. EMP_CPR_glb(1) == 0) THEN
                        DO k = 1,PD_Param_ID
                        ECOM_coef (k) = Xmatrix(9+k,1)
                        END DO

                        IF (NPARAM_glb /= PD_Param_ID + 3) THEN
                        PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
                        PRINT*,           'NPARAM_glb  =', NPARAM_glb
                        PRINT*,           'PD_Param_ID =', PD_Param_ID+3
                        PRINT*,'PROGRAM STOP AT m_orbdet.f03'
                        STOP
                        END IF

                END IF
        END IF


END IF

IF (ECOM_param_glb /= 0) THEN
   ECOM_accel_aposteriori = ECOM_accel_glb + ECOM_coef

!print*,'ECOM_coef=',ECOM_coef
!print*,'ECOM_accel_glb=',ECOM_accel_glb
!print*,'ECOM_accel_aposteriori=',ECOM_accel_aposteriori

! ----------------------------------------------------------------------
! SRP parameters
   ECOM_0_coef = ECOM_accel_aposteriori
!fname_id = PRN
   CALL doy2str(DOYSTR)
   fname_id = DOYSTR
   IF (ECOM_param_glb == 1) THEN
   fname = 'ECOM1_srp.in'
   param_id = 'ECOM1'
   write (param_value, *) ECOM_0_coef
   Call write_prmfile (fname, fname_id, param_id, param_value)
   END IF

   IF (ECOM_param_glb == 2) THEN
   fname = 'ECOM2_srp.in'
   param_id = 'ECOM2'
   write (param_value, *) ECOM_0_coef
   Call write_prmfile (fname, fname_id, param_id, param_value)
   END IF

   IF (ECOM_param_glb == 3) THEN
   fname = 'SBOXW_srp.in'
   param_id = 'SBOXW'
   write (param_value, *) ECOM_0_coef
   Call write_prmfile (fname, fname_id, param_id, param_value)
   END IF

END IF

! End of SRP model
! **********************************************************************


print *,"DELTA_V_apriori_glb", DELTA_V_apriori_glb  

! ----------------------------------------------------------------------
! Pseudo-stochastic pulses 
! ----------------------------------------------------------------------
IF (PULSE_param == 1) Then
	DO i_pulse = 1 , 3
		delta_v_corr(i_pulse) = Xmatrix(6 + NPARAM_EMP_ECOM_glb + i_pulse ,1)
	END DO
	DELTA_V_aposteriori_glb = DELTA_V_apriori_glb + delta_v_corr
END IF
! ----------------------------------------------------------------------

print *,"DELTA_V_aposteriori_glb", DELTA_V_aposteriori_glb

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write the estimated parameters in the input files
! ----------------------------------------------------------------------
write (fname_id, *) i
! ----------------------------------------------------------------------
! Initial state vector
If (1<0) Then
! SVEC_Zo_ESTIM transformation to ITRF

! or Reference_frame set to ICRF
fname = EQMfname
param_id = 'Reference_frame'
param_value = 'ICRF'
Call write_prmfile (fname, fname_id, param_id, param_value)
param_id = 'state_vector'
write (param_value, *) SVEC_Zo_ESTIM
Call write_prmfile (fname, fname_id, param_id, param_value)

fname = VEQfname
param_id = 'Reference_frame'
param_value = 'ICRF'
Call write_prmfile (fname, fname_id, param_id, param_value)
param_id = 'state_vector'
write (param_value, *) SVEC_Zo_ESTIM
Call write_prmfile (fname, fname_id, param_id, param_value)

End If
! ----------------------------------------------------------------------

End Do
! End of Orbit estimator iterations
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Integration :: Orbit final solution (after parameter estimation)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Determination and Prediction mode
IF (POD_MODE_glb == 2) THEN

! Orbit Determination number of epochs (without predicted orbit epochs)
sz1 = size(orb_icrf, DIM = 1)
sz2 = size(orb_icrf, DIM = 2)
Nepochs_estim = sz1

! ----------------------------------------------------------------------
! Rewrite Configuration files :: Add the Orbit Prediction arc length
! POD mode 4 : Orbit Determination and Orbit Prediction
! ----------------------------------------------------------------------
! Copy Configuration files 
fname_id = 'pred'
CALL write_prmfile2 (EQMfname, fname_id, EQMfname_pred)
CALL write_prmfile2 (VEQfname, fname_id, VEQfname_pred)

! Orbit arc length (estimated part) :: "orbarc" via module mdl_param.f03
!Call prm_main (EQMfname)

! Orbit arc length: Estimation + Prediction arc (in seconds)
orbarc_sum = orbarc + ORBPRED_ARC_glb
!print *,"POD_MODE_glb, EQMfname_pred ", POD_MODE_glb, EQMfname_pred
!print *,"orbit arc lengths           ", orbarc_sum, orbarc, ORBPRED_ARC_glb

param_id = 'Orbit_arc_length'
write (param_value, *) orbarc_sum
Call write_prmfile (EQMfname_pred, fname_id, param_id, param_value) 
Call write_prmfile (VEQfname_pred, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Integration 
! Numerical Integration: Variational Equations
VEQmode = 1
Call orbinteg (VEQfname_pred, VEQmode, orb0, veqSmatrix, veqPmatrix)
! Numerical Integration: Equation of Motion
VEQmode = 0
Call orbinteg (EQMfname_pred, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------

! Orbit estimated part without predicted part
ALLOCATE (orb_icrf_estim(Nepochs_estim,sz2), STAT = AllocateStatus)
if (AllocateStatus .ne. 0) then
        write(mesg, *) "Not enough memory - failed to allocate orb_icrf_estim, dimension=(", &
                Nepochs_estim, ",", sz2, ")"
        call report('FATAL', pgrm_name, 'orbdet', mesg, 'src/m_orbdet.f03', 1)
end if
orb_icrf_estim = orb_icrf(1:Nepochs_estim,1:sz2)

! Orbit residuals; statistics ! ICRF
Call statdelta(pseudobs_ICRF, orb_icrf_estim, dorb_icrf, RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr)
sz1 = size(dorb_icrf, DIM = 1)
sz2 = size(dorb_icrf, DIM = 2)
ALLOCATE (Vres(sz1,5), STAT = AllocateStatus)
if (AllocateStatus .ne. 0) then
        write(mesg, *) "Not enough memory - failed to allocate Vres, dimension=(", &
                sz1, ",5)"
        call report('FATAL', pgrm_name, 'orbdet', mesg, 'src/m_orbdet.f03', 1)
end if
Vres = dorb_icrf(1:sz1,1:5)
Vrms  = RMSdsr(1:3) 
!print *,"Orbit residuals opt (ICRF) RMS(XYZ)", RMSdsr(1:3)
!print *,"Orbit residuals: ICRF in XYZ" 
!WRITE (*,FMT='(A17, A4, 3F14.4)') "RMS-XYZ ICRF FIT", PRN, RMSdsr(1:3)

! Orbit residuals in orbital frame; statistics ! ICRF
!CALL statorbit (pseudobs_ICRF, orb_icrf_estim, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ, stat_RTN, stat_Kepler)
!print *,"Orbit residuals: ICRF in orbital frame" 
!WRITE (*,FMT='(A17, A4, 3F14.4)') "RMS-RTN ICRF FIT", PRN, stat_RTN(1, 1:3)
ELSE 

! ----------------------------------------------------------------------
! Orbit Integration 
VEQmode = 0
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit residuals; statistics ! ICRF
Call statdelta(pseudobs_ICRF, orb_icrf, dorb_icrf, RMSdsr, Sigmadsr, MEANdsr, MINdsr, MAXdsr)
sz1 = size(dorb_icrf, DIM = 1)
sz2 = size(dorb_icrf, DIM = 2)
ALLOCATE (Vres(sz1,5), STAT = AllocateStatus)
if (AllocateStatus .ne. 0) then
        write(mesg, *) "Not enough memory - failed to allocate Vres, dimension=(", &
                sz1, ",5)"
        call report('FATAL', pgrm_name, 'orbdet', mesg, 'src/m_orbdet.f03', 1)
end if
Vres = dorb_icrf(1:sz1,1:5)
Vrms  = RMSdsr(1:3)
!print *,"Orbit residuals opt (ICRF) RMS(XYZ)", RMSdsr(1:3)
!print *,"Orbit residuals: ICRF in XYZ" 
!WRITE (*,FMT='(A17, A4, 3F14.4)') "RMS-XYZ ICRF FIT", PRN, RMSdsr(1:3)

! Orbit residuals in orbital frame; statistics ! ICRF
!CALL statorbit (pseudobs_ICRF, orb_icrf, dorb_icrf, dorb_RTN, dorb_Kepler, stat_XYZ, stat_RTN, stat_Kepler)
!print *,"Orbit residuals: ICRF in orbital frame" 
!WRITE (*,FMT='(A17, A4, 3F14.4)') "RMS-RTN ICRF FIT", PRN, stat_RTN(1, 1:3)
! ----------------------------------------------------------------------

END IF

End If
! End Of Orbit Determination & Parameter estimation
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! POD modes :: 3, 4
! 3. Orbit Integration (Equation of Motion only)
! 4. Orbit Integration and Partials (Equation of Motion and Variational Equations)
! ----------------------------------------------------------------------
If (ESTmode == 0) then
! ----------------------------------------------------------------------
! POD mode 3
! ----------------------------------------------------------------------
IF      (VEQ_integration_glb == 0) THEN
! Numerical Integration: Equation of Motion
VEQmode = 0
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! POD mode 4
! ----------------------------------------------------------------------
ELSE IF (VEQ_integration_glb == 1) THEN
! Numerical Integration: Variational Equations
VEQmode = 1
Call orbinteg (VEQfname, VEQmode, orb0, veqSmatrix, veqPmatrix)
! Numerical Integration: Equation of Motion
VEQmode = 0
Call orbinteg (EQMfname, VEQmode, orb_icrf, veq0, veq1)
! ----------------------------------------------------------------------
END IF
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Control of orbit matrices dimensions in case that orbit integrator step
! has been reduced due to eclipse seasons 
! ----------------------------------------------------------------------
!CALL eclipse_integstep (EQMfname, VEQfname, mjd, r_sat, v_sat, integstep_flag, integstep_initial, integstep_reduced)
IF (integstep_flag) THEN
	integstep_rate = INT(integstep_initial/integstep_reduced) 
ELSE
	integstep_rate = 0
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit propagation backwards in order to propagate orbits and partials at epochs prior the ICs (initial conditions)
! ----------------------------------------------------------------------
IF (orbit_backwards_arc_cfg > 0.0D0) THEN
	orbintegr_back_flag = 1
ELSE
	orbintegr_back_flag = 0
END IF 

IF (orbintegr_back_flag > 0) THEN

! Copy Configuration files 
fname_id = 'back'
CALL write_prmfile2 (EQMfname, fname_id, EQMfname_back)
CALL write_prmfile2 (VEQfname, fname_id, VEQfname_back)

! Set negative sign to numerical integration step
IF (integstep_flag) THEN
	integstep_orb = integstep_reduced 
ELSE
	!integstep_orb = integstep_initial
	integstep_orb = integstep
END IF
integstep_orbback = -1.0D0 * integstep_orb
!integstep_orbback = -1.0D0 * integstep

param_id = 'integrator_step'
write (param_value, *)  integstep_orbback
Call write_prmfile (EQMfname_back, fname_id, param_id, param_value)
Call write_prmfile (VEQfname_back, fname_id, param_id, param_value)

! Set orbit arc length for backwards orbit propagation
! Orbit arc length (in seconds)
orbarc_back = orbit_backwards_arc_cfg * 3600.D0
param_id = 'Orbit_arc_length'
write (param_value, *) orbarc_back
Call write_prmfile (EQMfname_back, fname_id, param_id, param_value)
Call write_prmfile (VEQfname_back, fname_id, param_id, param_value)

! Orbit integration backwards: Equation of Motion solution
VEQmode = 0
Call orbinteg (EQMfname_back, VEQmode, orb_back, veq0, veq1)
! Orbit integration backwards: Variational Equations solution
VEQmode = 1
Call orbinteg (VEQfname_back, VEQmode, orb0, veqSmatrix_back, veqPmatrix_back)

! Merge orbits and partials matrices
! Orbit matrix
CALL matrixreverse (orb_back)
CALL matrixmerge   (orb_back, orb_icrf)
! veqSmatrix matrix
CALL matrixreverse (veqSmatrix_back)
CALL matrixmerge   (veqSmatrix_back, veqSmatrix)
! veqPmatrix matrix
CALL matrixreverse (veqPmatrix_back)
CALL matrixmerge   (veqPmatrix_back, veqPmatrix)

END IF  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Resize matrices; in case of no changes of the initial orbit integrator step
! the output matrix is a copy of the input one
! ----------------------------------------------------------------------
! Orbit Matrix 
CALL orbresize (orb_icrf, integstep_rate, orb_icrf_final) 
! Partials matrices
IF (ESTmode > 1 .OR. VEQ_integration_glb == 1) THEN
CALL orbresize (veqSmatrix, integstep_rate, veqSmatrix_final) 
CALL orbresize (veqPmatrix, integstep_rate, veqPmatrix_final) 
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit transformation to terrestrial frame: ICRF to ITRF
! ----------------------------------------------------------------------
! Time System according to global variable TIME_Scale (Module mdl_param.f03)
time_sys = TIME_SCALE
!CALL orbC2T (orb_icrf, time_sys, orb_itrf)
CALL orbC2T (orb_icrf_final, time_sys, orb_itrf_final)
! ----------------------------------------------------------------------

 100  if (allocated(ECOM_0_coef)) deallocate(ECOM_0_coef, stat=DeallocateStatus)
      END SUBROUTINE

End

