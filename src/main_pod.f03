      program main_pod


! ----------------------------------------------------------------------
! Program:	main_pod.f03
! ----------------------------------------------------------------------
! Purpose:
!  Precise Orbit Determination (POD) of GNSS satellites 
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
! - Dr. Thomas Papanikolaou, 29 March 2019
! 	POD version upgrade to a multi-GNSS multi-satellite POD version 
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_config
      USE mdl_param
      USE m_read_leapsec
      USE m_pod_gnss
      USE m_writeorbit_multi
      USE m_writearray
      USE m_writearray2
      USE m_writeorbit
	  USE m_write_orb2sp3
	  USE m_clock_read
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: EQMfname, VEQfname, PODfname, ORBMODfname				
      CHARACTER (LEN=100) :: EQMfname_initial, VEQfname_initial				
      INTEGER (KIND = prec_int2) :: ios_line, ios_key, ios_data
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbits_ics_icrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Vres  
      REAL (KIND = prec_d), DIMENSION(3) :: Vrms 	    
	  !REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT
! ----------------------------------------------------------------------
      CHARACTER (LEN=2) :: GNSS_id
	  INTEGER (KIND = prec_int2) :: ORB_mode
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Nsat, isat
	  INTEGER (KIND = prec_int8) :: iepoch, iparam
	  INTEGER (KIND = prec_int8) :: i
	  INTEGER (KIND = prec_int8) :: sz1, sz2, Nepochs, N2_orb, N2_veqSmatrix, N2_veqPmatrix, N2sum  
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_icrf  
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_itrf  
	  CHARACTER (LEN=3), ALLOCATABLE :: PRNmatrix(:)
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus  
	  CHARACTER (LEN=3) :: PRN_isat
	  INTEGER :: ios
      CHARACTER (LEN=100) :: orbits_fname, orbits_partials_fname				
      CHARACTER (LEN=100) :: fname_write				
      CHARACTER (LEN=100) :: filename				
! ----------------------------------------------------------------------
      CHARACTER (LEN=300) :: fname_sp3, ORBpseudobs_fname, ORBEXT_fname				
	  INTEGER :: year, month, day
	  INTEGER :: Iyear, Imonth, Iday
      REAL (KIND = prec_d) :: Sec_00 	    
! ----------------------------------------------------------------------
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
      REAL (KIND = prec_d) :: Zo(6) 
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: ORB2sp3_fname				
      INTEGER (KIND = prec_int2) :: sat_vel	  	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resR  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resT  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resN  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: mjd
      INTEGER (KIND = prec_int8) :: GPS_week, GPSweek_mod1024
      REAL (KIND = prec_d) :: GPS_wsec, GPS_day
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Ncommon  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_RTN, dorb_Kepler
! ----------------------------------------------------------------------
      !INTEGER (KIND = prec_int2) :: POD_MODE_glb	  	  
      !REAL (KIND = prec_d) :: ORBPRED_ARC_glb
      REAL (KIND = prec_d) :: orbarc_sum
      REAL (KIND = prec_d) :: orb_est_arc
      INTEGER (KIND = prec_int2) :: IC_MODE	  	  
      CHARACTER (LEN=500) :: IC_REF				
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbdiff2
! ----------------------------------------------------------------------
      LOGICAL :: pod_config_exists
	  CHARACTER (LEN=100) :: pgm_name
! ----------------------------------------------------------------------
      CHARACTER (len=300) :: str
      INTEGER (KIND = prec_int2) :: j
      CHARACTER (len=9) :: POD_version
      REAL (KIND = prec_q), DIMENSION(:,:,:), ALLOCATABLE :: CLKmatrix 
      CHARACTER (LEN=300) :: CLKfname
      INTEGER (KIND = prec_int2) :: CLKformat
! ----------------------------------------------------------------------

NPARAM_glb = 0

! CPU Times
CALL cpu_time (CPU_t0)

! ----------------------------------------------------------------------
! POD Version:
POD_version = 'v.1.0.1'
! ----------------------------------------------------------------------

! init globals
CALL globals_init()

! ----------------------------------------------------------------------
! POD major configuration file
! ----------------------------------------------------------------------
! Default master POD configureation file
PODfname = 'POD.in'

! Read command line to see if non default master configration file given
CALL read_cmdline

! Check if non-default config file given on the command line
If ( trim(POD_fname_cfg) .ne. 'DEFAULT' ) then
	PODfname = trim(POD_fname_cfg)
End If

! Check for existance of POD config file
pod_config_exists = .true.
INQUIRE(FILE=PODfname, EXIST=pod_config_exists)

If ( .not. pod_config_exists) then
	call get_command_argument( 0, pgm_name )
    write(*,'(3a)') 'No Default config file found (POD.in)  - Type: ',trim(pgm_name),' --help'
    write(*,'(3a)') 'If using a non-default config.filename - Type: ',trim(pgm_name),' -c config.filename'
	STOP
End If

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------  ! 999999999999999999999999999
! ----------------------------------------------------------------------
! Temporary:: manual configuration 
! ----------------------------------------------------------------------
! Beidou case:
! 1. BDSorbtype = 'IGSO'
! 2. BDSorbtype = 'MEO'
BDSorbtype_glb = 'MEO'

! Empirical Forces reference frame:
! 1. Orbital Frame
! 2. Body-fixed frame
Frame_EmpiricalForces_glb = 1

!print *,"Frame_EmpiricalForces_glb ", Frame_EmpiricalForces_glb
!print *,"BDSorbtype_glb            ", BDSorbtype_glb
!print *,"              "
! ----------------------------------------------------------------------
! Temporary :: End of input POD configuration
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------  ! 999999999999999999999999999

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Start :: Read major configuration file POD.in
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! POD Tool mode:
! ----------------------------------------------------------------------
! 1. Orbit Determination (pseudo-observations; orbit fitting)
! 2. Orbit Determination and Prediction
! 3. Orbit Integration (Equation of Motion only)
! 4. Orbit Integration and Partials (Equation of Motion and Variational Equations)
! ----------------------------------------------------------------------
param_id = 'POD_MODE_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) POD_MODE_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Initial Conditions input mode
! ----------------------------------------------------------------------
! 1. Input a-priori orbit in sp3 format (applied as pseudo-observations)
! 2. Input file with Initial Conditions (State Vector and Parameters at initial epoch per satellite) 
! ----------------------------------------------------------------------
param_id = 'IC_input'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) IC_MODE_cfg 

! Initial Conditions reference frame
param_id = 'IC_refsys'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) IC_REF_cfg

! Initial Conditions file name
param_id = 'IC_filename_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) IC_filename_cfg
! ----------------------------------------------------------------------
!print *,"IC_MODE_cfg, IC_REF_cfg", IC_MODE_cfg, IC_REF_cfg
!print *,"IC_filename_cfg", IC_filename_cfg

! ----------------------------------------------------------------------
! Configuration files of Orbit modelling (2 Basic initial files):
! ----------------------------------------------------------------------
! Equation of Motion
param_id = 'EQM_fname_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) EQM_fname_cfg 

! Variational Equations
param_id = 'VEQ_fname_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) VEQ_fname_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! POD_MODE Cases: 1 or 2  
! ----------------------------------------------------------------------
! A-priori orbit sp3 as pseudo-observations :: sp3 file name
! ----------------------------------------------------------------------
param_id = 'pseudobs_orbit_filename_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) pseudobs_orbit_filename_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! POD_MODE Cases: ALL  (IF orbit_external_opt .NE. 0 see EQM.in configuration file)
! ----------------------------------------------------------------------
! External Orbit Comparison :: sp3 file name
! ----------------------------------------------------------------------
param_id = 'ext_orbit_filename_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) ext_orbit_filename_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Determination arc length
! ----------------------------------------------------------------------
param_id = 'orbit_determination_arc_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) orbit_determination_arc_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! POD_MODE Cases: 2  
! ----------------------------------------------------------------------
! Orbit Prediction arc length (Seconds)
! ----------------------------------------------------------------------
param_id = 'orbit_prediction_arc_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) orbit_prediction_arc_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Propagation backwards arc length
! ----------------------------------------------------------------------
param_id = 'orbit_backwards_arc_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) orbit_backwards_arc_cfg 
! ----------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Earth Orientation Parameters (EOP)
! ---------------------------------------------------------------------------
! EOP data solution options:
! 1. IERS C04 										: EOP_sol=1
! 2. IERS RS/PC Daily 								: EOP_sol=2
! 3. IGS ultra-rapid ERP + IERS RS/PC Daily (dX,dY)	: EOP_sol=3
param_id = 'EOP_solution_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) EOP_solution_cfg 

! EOP filename by IERS EOP :: Solutions 1 and 2
param_id = 'EOP_fname_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) EOP_fname_cfg 

! ERP filename (Earth Rotation Parameters by IGS) :: Solution 3
param_id = 'ERP_fname_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) ERP_fname_cfg 

! EOP data interpolation number of points	  
param_id = 'EOP_Nint_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) EOP_Nint_cfg 
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! IAU Precession-Nutation model:
! ---------------------------------------------------------------------------
! 1. IAU2000A:		iau_pn_model = 2000
! 2. IAU2006/2000A:	iau_pn_model = 2006
param_id = 'iau_model_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) iau_model_cfg 
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Orbit Parameter Estimation
! ---------------------------------------------------------------------------
! Number of iterations
param_id = 'Estimator_Iterations_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) Estimator_Iterations_cfg 
! ---------------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write to sp3 orbit format: Option for write Satellite Velocity vector 
! ----------------------------------------------------------------------
! 0. sat_vel = 0 :: Do not write Velocity vector to sp3 orbit
! 1. sat_vel > 0 :: Write Velocity vector to sp3 orbit
param_id = 'sp3_velocity_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) sp3_velocity_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write partials of the velocity vector w.r.t. parameters into the orbits_partials output file 
! ----------------------------------------------------------------------
! 0. partials_velocity_cfg = 0 :: Do not write Velocity vector's partials elements
! 1. partials_velocity_cfg > 0 :: Write Velocity vector's partials elements
param_id = 'partials_velocity_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) partials_velocity_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Leap Second file name:
! ----------------------------------------------------------------------
param_id = 'leapsec_filename_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) leapsec_filename_cfg 
! ----------------------------------------------------------------------

! Read Satellite infromation from SINEX file
! ----------------------------------------------------------------------
param_id = 'satsinex_filename_cfg'
Call readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) satsinex_filename_cfg
PRINT*,'satsinex_filename =', satsinex_filename_cfg
!----------------------------------------------------------------------

! Flag of different models for A priori SRP value 
! ---------------------------------------------------------------------
param_id = 'Flag_BW_cfg'
Call readparam (PODfname, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) Flag_BW_cfg
PRINT*,'Flag_BW_cfg =', Flag_BW_cfg
! ----------------------------------------------------------------------
! Reference System of Variational Equations' Partials & Parameter Estimation 
! ----------------------------------------------------------------------
! 1. Celestial Reference System :: ICRS
! 2. Terrestrial Reference System :: ITRS
! ----------------------------------------------------------------------
param_id = 'VEQ_REFSYS_cfg'
CALL readparam (PODfname, param_id, param_value)
READ ( param_value, FMT ='(A4)', IOSTAT=ios_key ) VEQ_REFSYS_cfg 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! End :: Read major configuration file POD.in
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! Read command line to oerwrite POD.in configfile options
CALL read_cmdline
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Major configuration parameters via "read_cmdline routine", "POD.in configuration file" and "Module mdl_config.f03" 
! ----------------------------------------------------------------------
! Debug
if (1<0) then
print*, 'Master configuration parameters:'
print*, 'POD_MODE_cfg: '               ,POD_MODE_cfg
print*, 'EQM_fname_cfg: '              ,trim(EQM_fname_cfg)
print*, 'VEQ_fname_cfg: '              ,trim(VEQ_fname_cfg)
print*, 'pseudobs_orbit_filename_cfg: ',trim(pseudobs_orbit_filename_cfg) 
print*, 'ext_orbit_filename_cfg: '     ,trim(ext_orbit_filename_cfg)
print*, 'orbit_determination_arc_cfg: ',orbit_determination_arc_cfg
print*, 'orbit_prediction_arc_cfg: '   ,orbit_prediction_arc_cfg
print*, 'EOP_solution_cfg: '           ,EOP_solution_cfg
print*, 'EOP_fname_cfg: '              ,trim(EOP_fname_cfg)
print*, 'ERP_fname_cfg: '              ,trim(ERP_fname_cfg)
print*, 'EOP_Nint_cfg: '               ,EOP_Nint_cfg
print*, 'iau_model_cfg: '              ,iau_model_cfg
print*, 'Estimator_Iterations_cfg: '   ,Estimator_Iterations_cfg 
print*, 'sp3_velocity_cfg: '           ,sp3_velocity_cfg
print*, 'leapsec_filename_cfg: '        ,leapsec_filename_cfg
end if
! ----------------------------------------------------------------------

POD_MODE_glb      = POD_MODE_cfg
orb_est_arc       = orbit_determination_arc_cfg * 3600.D0
ORBPRED_ARC_glb   = orbit_prediction_arc_cfg * 3600.D0

! Debug
if (1<0) then
print *, "POD_MODE_glb", POD_MODE_glb
print *, "orb_est_arc", orb_est_arc
print *, "ORBPRED_ARC_glb", ORBPRED_ARC_glb
Print *," "
end if

! POD mode
If      (POD_MODE_glb == 1) then
Print *,"POD Tool mode: 1 :: Orbit Determination"
Else IF (POD_MODE_glb == 2) then 
Print *,"POD Tool mode: 2 :: Orbit Determination and Prediction"
Else IF (POD_MODE_glb == 3) then 
Print *,"POD Tool mode: 3 :: Orbit Integration"
Else IF (POD_MODE_glb == 4) then 
Print *,"POD Tool mode: 4 :: Orbit Integration and Partials"
End IF
Print *," "
! IC mode
If      (IC_MODE_cfg == 1) then
Print *,"Initial Conditions mode: 1 :: A-priori orbit input file: ", TRIM(pseudobs_orbit_filename_cfg)
Else IF (IC_MODE_cfg == 2) then 
Print *,"Initial Conditions mode: 2 :: Initial Conditions input file: ", TRIM(IC_filename_cfg)
END IF
Print *," "

! ----------------------------------------------------------------------
! Read Leap Second File
CALL read_leapsec(leapsec_filename_cfg)
!Print*,'NDAT,IDAT,DATS: ',NDAT,IDAT,DATS

! ----------------------------------------------------------------------
! Form (rewrite) the two orbit integration configuration files for 
! Equation of Motion and Variational Equations: EQM.in and VEQ.in 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Copy Initial Configuration files 
! ----------------------------------------------------------------------
fname_id = '0'
CALL write_prmfile2 (EQM_fname_cfg, fname_id, EQMfname)
CALL write_prmfile2 (VEQ_fname_cfg, fname_id, VEQfname)
! ----------------------------------------------------------------------
!print *,"EQMfname ", EQMfname
!print *,"VEQfname ", VEQfname
!print *,"              "

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Start :: Rewrite configuration files
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set POD_MODE
! ----------------------------------------------------------------------
fname_id = '1'
IF (POD_MODE_glb == 1 .OR. POD_MODE_glb == 2) THEN
param_id = 'VEQ_integration'
param_value = '1'
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'Estimator_procedure'
param_value = '1'
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

ELSE IF (POD_MODE_glb == 3) THEN
param_id = 'VEQ_integration'
param_value = '0'
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'Estimator_procedure'
param_value = '0'
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

ELSE IF (POD_MODE_glb == 4) THEN
param_id = 'VEQ_integration'
param_value = '1'
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'Estimator_procedure'
param_value = '0'
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set orbit estimation arc length
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'Orbit_arc_length'
write (param_value, *) orb_est_arc
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set a-priori orbit (pseudo-observations; orbit comparison)
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'pseudobs_filename'
param_value = pseudobs_orbit_filename_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'orbit_filename'
param_value = ext_orbit_filename_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set Earth Orientation Modelling 
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'EOP_data_sol'
write (param_value, *) EOP_solution_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'EOP_filename'
write (param_value, *) EOP_fname_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'ERP_filename'
write (param_value, *) ERP_fname_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'EOP_interpolation_points'
write (param_value, *) EOP_Nint_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'iau_pn_model'
write (param_value, *) iau_model_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set Orbit Parameter Estimator number of iterations 
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'Estimator_Iterations'
write (param_value, *) Estimator_Iterations_cfg
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! End :: Rewrite configuration files
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! POD of the GNSS satellite constellations
! ----------------------------------------------------------------------
CALL pod_gnss (EQMfname, VEQfname, PRNmatrix, orbits_partials_icrf, orbits_partials_itrf, &
               orbits_ics_icrf,orbit_resR, orbit_resT, orbit_resN, orbdiff2)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Clocks read from input file for passing to the write out 
! ----------------------------------------------------------------------
IF (IC_MODE_cfg == 1) THEN 
CLKformat = 1
CLKfname  = pseudobs_orbit_filename_cfg
ELSE
CLKformat = 0
CLKfname = ''
END IF 
CALL clock_read (CLKfname,CLKformat, PRNmatrix, CLKmatrix)
! ---------------------------------------------------------------------- 

! ----------------------------------------------------------------------
! Output filenames prefix
! ----------------------------------------------------------------------
!mjd = orbits_partials_icrf(1,1,1)
!mjd = orbits_partials_itrf(2,1,1)
mjd = orbits_ics_icrf(1,1) 
CALL time_GPSweek  (mjd, GPS_week, GPS_wsec, GPSweek_mod1024)
!CALL time_GPSweek2 (mjd, GPS_week, GPS_wsec, GPSweek_mod1024, GPS_day)
GPS_day = ( GPS_wsec/86400.0D0 )  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write satellite orbits and partial derivatives to one .orb output file (internal format)
! ----------------------------------------------------------------------
!orbits_partials_fname = 'orbits_partials_icrf.orb'
write (orbits_partials_fname, FMT='(A3,I4,I1,A20)') 'gag', (GPS_week), INT(GPS_day) ,'_orbits_partials.out'
!CALL writeorbit_multi (orbits_partials_icrf, PRNmatrix, orbits_partials_fname)
CALL writeorbit_multi (orbits_partials_icrf, orbits_partials_itrf, orbits_ics_icrf, PRNmatrix, & 
						orbits_partials_fname, EQMfname, VEQfname, POD_version)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write satellite orbits to sp3 format
! ----------------------------------------------------------------------
! Orbit sp3 filename
write (ORB2sp3_fname, FMT='(A3,I4,I1,A4)') 'gag', (GPS_week), INT(GPS_day) ,'.sp3'
sat_vel = sp3_velocity_cfg
! ICRF
!CALL write_orb2sp3 (orbits_partials_icrf, PRNmatrix, ORB2sp3_fname, sat_vel, CLKmatrix)
! ITRF
CALL write_orb2sp3 (orbits_partials_itrf, PRNmatrix, ORB2sp3_fname, sat_vel, CLKmatrix)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Extract residual filename tag from input .sp3 filename
! ----------------------------------------------------------------------
str = trim(adjustl(pseudobs_orbit_filename_cfg))
i = index(str, '.sp3')
j = len(str(1:i-1))

! ----------------------------------------------------------------------
! Write Orbit residuals
! ----------------------------------------------------------------------
! Radial
write (filename, FMT='(A3,I4,I1,a1,a,A16)') 'gag', (GPS_week), INT(GPS_day), '_', str(1:j), '_orbitstat_R.out'
Call writearray (orbit_resR, filename)
! Transverse
write (filename, FMT='(A3,I4,I1,a1,a,A16)') 'gag', (GPS_week), INT(GPS_day), '_', str(1:j) ,'_orbitstat_T.out'
Call writearray (orbit_resT, filename)
! Normal
write (filename, FMT='(A3,I4,I1,a1,a,A16)') 'gag', (GPS_week), INT(GPS_day), '_', str(1:j) ,'_orbitstat_N.out'
Call writearray (orbit_resN, filename)
! ----------------------------------------------------------------------
! Write combined orbit residuals file (RTN)
write (filename, FMT='(A3,I4,I1,a1,a,A16)') 'gag', (GPS_week), INT(GPS_day), '_', str(1:j) ,'_orbdiff_rtn.out'
Call writearray2 (orbdiff2, filename)

CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0

End Program

