MODULE m_writeorbit_multi


! ----------------------------------------------------------------------
! MODULE: m_writeorbit_multi.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for write orbit array data to output (ascii) files 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	21 March 2019
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			

	  
Contains


SUBROUTINE writeorbit_multi (orbitsmatrix_crf,orbitsmatrix_trf,orbits_ics_icrf,PRN_array,filename,EQMfname,VEQfname,POD_version)

! ----------------------------------------------------------------------
! SUBROUTINE: writeorbit_multi 
! ----------------------------------------------------------------------
! Purpose:
!  Write orbit and partial derivatives matrices to an output ascii file
!
!   writeorbit.f03 subroutine has been modified in order to write the 
!	estimated orbits and partial derivatives (solution of the Variational Equations) 
!   to output file (ascii) based on an internal adopted orbit format:
!   {MJD Sec_00h r(XYZ) v(XYZ) State_Transition_Matrix Sensitivity_matrix} 
! ----------------------------------------------------------------------
! Input arguments:
! - wrtArray:       Input allocatable array
! - filename:       Orbits & Partials file name to be used for writing out the orbits/partials matrices 
! - EQMfname: 		Input configuration file name for the orbit integration of the Equation of Motion  
! - VEQfname: 		Input configuration file name for the orbit integration of the Variational Equations
!
! Output arguments:
!
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	21 March 2019
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_config
      USE mdl_param
      USE m_read_satsnx
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:,:), ALLOCATABLE :: orbitsmatrix_crf
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:,:), ALLOCATABLE :: orbitsmatrix_trf
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:), ALLOCATABLE :: orbits_ics_icrf
      CHARACTER (LEN=3), ALLOCATABLE :: PRN_array(:)
      CHARACTER (LEN=100), INTENT(IN) :: filename
      CHARACTER (LEN=100), INTENT(IN) :: EQMfname, VEQfname
      CHARACTER (len=9), INTENT(IN) :: POD_version	  
! OUT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, i_write
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios, ios_ith, ios_key
      INTEGER (KIND = prec_int8) :: sz1, sz2, sz3
      INTEGER (KIND = prec_int2) :: wrt_opt
      INTEGER (KIND = prec_int2) :: FMT_opt
! ----------------------------------------------------------------------
      CHARACTER (LEN=1) :: RealT
      INTEGER (KIND = prec_int2) :: RealW, RealD
      CHARACTER (LEN=70) :: fmt_wrt, fmt_wrt0, fmt_sz2
      REAL (KIND = prec_q) :: wrtArrayLN 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: Nepochs, Nparam, Nsat, Norbits_ics_icrf 
      INTEGER (KIND = prec_int8) :: i_epoch, i_sat
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
      CHARACTER (LEN=100) :: gravity_model_filename, iau_pn_model, DE_fname_data, ocean_tides_model_file
      CHARACTER (LEN=100) :: orbit_num_integrator, EOP_sol, EOP_data, SE_Tides, Pole_Tide 
      CHARACTER (LEN=5) :: EPH_name
      INTEGER (KIND = prec_int4) :: ln
      INTEGER (KIND = prec_int2) :: Ftides
! ----------------------------------------------------------------------
      CHARACTER (LEN=3) :: PRN_isat
      !INTEGER (KIND = prec_int4) :: IY, IM, ID
      INTEGER Iyear, Imonth, Iday, J_flag
      DOUBLE PRECISION FD  
      REAL (KIND = prec_d) :: Sec_00, mjd, mjd_1, jd0
      !INTEGER (KIND = prec_int4) :: DOY
      CHARACTER (LEN=10) :: srp_model
! ----------------------------------------------------------------------
    CHARACTER(LEN=8)  :: date_mach
    CHARACTER(LEN=10) :: time_mach
    CHARACTER(LEN=5)  :: zone_mach
    INTEGER :: time_values(8)
    REAL (KIND = prec_q) :: kepler_ic(9),r_ic(3),v_ic(3)
! ----------------------------------------------------------------------


UNIT_IN = 7  												

! ----------------------------------------------------------------------
! Orbit arrays dimensions
sz1 = SIZE (orbitsmatrix_crf,DIM=1)
sz2 = SIZE (orbitsmatrix_crf,DIM=2)
sz3 = SIZE (orbitsmatrix_crf,DIM=3)

Nepochs = sz1
Nparam  = sz2
Nsat    = sz3
   
! PRN
Nsat = SIZE (PRN_array,DIM=1)

! orbits_ics_icrf
Norbits_ics_icrf = SIZE (orbits_ics_icrf,DIM=1)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Format definition
! ----------------------------------------------------------------------
! Orbit format: {MJD Sec_00h r(XYZ) v(XYZ)} 
!fmt_wrt = '(F25.15,F25.9,3F25.3,3F25.7)'
! Orbit-VEQ format: {PRN MJD Sec_00h r(XYZ) v(XYZ) VEQ-Z VEQ-P} 
!fmt_wrt = '(5A3,F25.12,F25.9,3F25.4,3F25.9, F25)'
fmt_wrt = '(A3,A1,F25.12,F25.9,3F25.4,3F25.9, A)'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Open file
      OPEN (UNIT=UNIT_IN,FILE=filename,ACTION="WRITE",POSITION="REWIND", IOSTAT=ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", filename
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write file header information
! ----------------------------------------------------------------------
WRITE (UNIT=UNIT_IN,FMT='(A,17x,A)',IOSTAT=ios_ith) '#INFO    POD Tool version :',TRIM(POD_version)
WRITE (UNIT=UNIT_IN,FMT='(A,17x,A)',IOSTAT=ios_ith) '#INFO    POD Tool output  :','Orbits & Partial Derivatives file'
CALL date_and_time(date_mach,time_mach,zone_mach,time_values)
WRITE (UNIT=UNIT_IN,FMT='(A,2x,I4,1x,I2,1x,I2,1x, I2,1x,I2,1x,I2,1x)',IOSTAT=ios_ith) &
       &'#INFO    File creation date (y/m/d/h/m/s):', &
       & time_values(1),time_values(2),time_values(3), time_values(5),time_values(6),time_values(7)      
	
! POD Tool mode
IF (POD_MODE_cfg == 1) WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    POD Tool mode:                     ', & 
					   & '1. Orbit Determination (using pseudo-observations)'
IF (POD_MODE_cfg == 2) WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    POD Tool mode:                     ', & 
					   & '2. Orbit Determination and Prediction'
IF (POD_MODE_cfg == 3) WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    POD Tool mode:                     ', & 
					   & '3. Orbit Numerical Integration'
IF (POD_MODE_cfg == 4) WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    POD Tool mode:                     ', & 
					   & '4. Orbit & Partials numerical integration '
! Initial Conditions
IF (IC_MODE_cfg == 1) WRITE (UNIT=UNIT_IN,FMT='(A,A,A)',IOSTAT=ios_ith) '#INFO    Initial Conditions input mode:     ', & 
					   & '1. a-priori .sp3 orbit: ',TRIM(pseudobs_orbit_filename_cfg)
IF (IC_MODE_cfg == 2) WRITE (UNIT=UNIT_IN,FMT='(A,A,A)',IOSTAT=ios_ith) '#INFO    Initial Conditions input mode:     ', & 
					   & '2. Initial Conditions input file: ',TRIM(IC_filename_cfg) 

WRITE (UNIT=UNIT_IN,FMT='(a,I9,F25.10)' ,IOSTAT=ios_ith) '#INFO    Epoch initial conditions:          ', &
                                                         INT(orbits_ics_icrf(1,1)),orbits_ics_icrf(2,1) !INT(orbitsmatrix_crf(1,1,1)),orbitsmatrix_crf(1,2,1)
WRITE (UNIT=UNIT_IN,FMT='(a,I9,F25.10)' ,IOSTAT=ios_ith) '#INFO    Epoch Start:                       ', & 
                                                         INT(orbitsmatrix_crf(1,1,1)),orbitsmatrix_crf(1,2,1) 
WRITE (UNIT=UNIT_IN,FMT='(a,I9,F25.10)' ,IOSTAT=ios_ith) '#INFO    Epoch End:                         ', &
                                                         INT(orbitsmatrix_crf(Nepochs,1,1)),orbitsmatrix_crf(Nepochs,2,1) 
WRITE (UNIT=UNIT_IN,FMT='(a,I9)'        ,IOSTAT=ios_ith) '#INFO    Tabular interval (sec):            ', abs(INT(integstep))
WRITE (UNIT=UNIT_IN,FMT='(a,i5)'        ,IOSTAT=ios_ith) '#INFO    Number of Epochs:                  ',Nepochs

! Orbit arc length (in hours) 
WRITE (UNIT=UNIT_IN,FMT='(A,A,I3,A,I3,A,I3)',IOSTAT=ios_ith) '#INFO    Orbit arc length (hours):         ',& 
								  'Orbit Determination arc: ', INT(orbit_determination_arc_cfg), &
								  ' | Orbit Prediction arc: ', INT(orbit_prediction_arc_cfg), &
								  ' | Backwards orbit integration arc: ', INT(orbit_backwards_arc_cfg)

! Numerical Integration methods
IF (integmeth == 1) orbit_num_integrator = 'RKN7(6)8 Runge-Kutta-Nystrom 7th order method'
IF (integmeth == 2) orbit_num_integrator = 'Runge-Kutta 4th order'
IF (integmeth == 3) orbit_num_integrator = 'RK8(7)13 Runge-Kutta 8th order'
WRITE (UNIT=UNIT_IN,FMT='(2A)',IOSTAT=ios_ith)   '#INFO    Numerical Integration Method:     ', orbit_num_integrator 
WRITE (UNIT=UNIT_IN,FMT='(A,I7)',IOSTAT=ios_ith) '#INFO    Numerical Integration step (sec): ', abs(INT(integstep)) 

! Satellite Dynamics model
!WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO    Model information:           [TIME_SYS] [GRAV_MODEL]&
!                                           &[PN_MODEL] [EPH_MODEL] [ALBEDO_MODEL]'
WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    Time System:                       ', TIME_SCALE 
param_id = 'gravity_model_filename'
CALL readparam (EQM_fname_cfg, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) gravity_model_filename 
param_id = 'DE_fname_data'
CALL readparam (EQM_fname_cfg, param_id, param_value)
!READ ( TRIM(ADJUSTR(param_value)), FMT = * , IOSTAT=ios_key ) DE_fname_data
DE_fname_data =  TRIM(param_value) 
ln = LEN_TRIM (DE_fname_data)
WRITE (EPH_name, FMT = '(A2,A3)' , IOSTAT=ios_key ) 'DE', DE_fname_data(ln-2:ln)

! General orbit parameterization											
Call prm_main (EQMfname)

IF (FMOD_GRAV(3) == 1 .and. FMOD_TIDES(1) == 1) THEN
SE_Tides = 'IERS 2010'
ELSE
SE_Tides = '-'
END IF
IF (FMOD_GRAV(3) == 1 .and. FMOD_TIDES(3) == 1) THEN
param_id = 'ocean_tides_model_file'
CALL readparam (EQM_fname_cfg, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) ocean_tides_model_file 
ELSE
ocean_tides_model_file = '-'
END IF
IF (FMOD_GRAV(3) == 1 .and. FMOD_TIDES(4) == 1) THEN
Pole_Tide = 'IERS 2010'
ELSE
Pole_Tide = '-'
END IF
WRITE (UNIT=UNIT_IN,FMT='(A)',IOSTAT=ios_ith)   '#INFO    Satellite Dynamics Model::         ' 
WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    Gravity field model:               ', TRIM(gravity_model_filename)
WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    Planetary Ephemeris:               ', TRIM(EPH_name)
WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    Solid Earth Tides  :               ', TRIM(SE_Tides)
WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    Ocean Tides        :               ', TRIM(ocean_tides_model_file)
WRITE (UNIT=UNIT_IN,FMT='(A,A)',IOSTAT=ios_ith) '#INFO    Pole Tide          :               ', TRIM(Pole_Tide)

! Earth Orientation modelling
IF (EOP_solution_cfg == 1) THEN 
	EOP_sol = 'IERS C04'
	EOP_data = EOP_fname_cfg
ELSE IF (EOP_solution_cfg == 2) THEN
	EOP_sol = 'IERS RS/PC Daily'
	EOP_data = EOP_fname_cfg
ELSE IF (EOP_solution_cfg == 3) THEN
	EOP_sol = 'IGS ultra-rapid ERP'
	EOP_data = ERP_fname_cfg
END IF
param_id = 'iau_pn_model'
CALL readparam (EQM_fname_cfg, param_id, param_value)
READ ( param_value, FMT = * , IOSTAT=ios_key ) iau_pn_model 
WRITE (UNIT=UNIT_IN,FMT='(10A)'            ,IOSTAT=ios_ith) '#INFO    Earth Orientation modelling:       ', &
											& '[EOP solution: ',trim(EOP_sol),'] ', &
											& '[EOP data file: ',trim(EOP_data),'] ', &
											& '[Precession-Nutation model: IAU',TRIM(iau_pn_model),']'

WRITE (UNIT=UNIT_IN,FMT='(a,i3)'           ,IOSTAT=ios_ith) '#INFO    Number of Satellites:              ',Nsat
WRITE (UNIT=UNIT_IN,FMT='(a,i4)'           ,IOSTAT=ios_ith) '#INFO    Number of Parameters per satellite:',NPARAM_glb+6
WRITE (UNIT=UNIT_IN,FMT='(a,i4)'           ,IOSTAT=ios_ith) '#INFO    Number of Partials:                ',Nparam-8
WRITE (UNIT=UNIT_IN,FMT='(a,a4)'           ,IOSTAT=ios_ith) '#INFO    Partials Reference System:         ',VEQ_REFSYS_cfg
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO    Satellite ICS:                     '
DO i_sat = 1 , Nsat
   ! Read SINEX file for SVN number, Mass, ..  
   !mjd0   = INT(orbits_ics_icrf(1,i_sat))
   Sec_00 = orbits_ics_icrf(2,i_sat)
   mjd = orbits_ics_icrf(1,i_sat)
   jd0 = 2400000.5D0
   CALL iau_JD2CAL ( jd0, mjd, Iyear, Imonth, Iday, FD, J_flag )
   CALL iau_CAL2JD ( Iyear, 1, 1, jd0, mjd_1, J_flag )   
   DOY = INT(mjd) - (mjd_1-1) 
   PRN_isat = PRN_array(i_sat)
   CALL read_satsnx (satsinex_filename_cfg, Iyear, DOY, Sec_00, PRN_isat) 
   ! SRP model
	IF (ECOM_param_glb == 1) THEN
      srp_model = 'ECOM1'   
	ELSE IF (ECOM_param_glb == 2) THEN
      srp_model = 'ECOM2'   
        ELSE IF (ECOM_param_glb == 3) THEN
      srp_model = 'SBOXW'
	END IF
! IC INFO   
   WRITE (UNIT=UNIT_IN,FMT='(a,a,1x,a3,1x,a,1x,i3,1x,a,1x,a,1x,a,1x,F10.5,1x,a,1x,a,1x,a,1x,i3,1x,a)' ,IOSTAT=ios_ith) & 
          &'#IC_INFO ','PRN:',PRN_array(i_sat),'SVN:',SVNID,'BLK_TYP:',TRIM(BLKTYP),' MASS:',MASS, &
          &'SRP:', TRIM(srp_model), 'Nparam:', NPARAM_glb+6, '[IC PARAM LIST - X Y Z XD YD ZD DRAD .....]'
! IC 
   WRITE (UNIT=UNIT_IN,FMT='(a,a3,1x,I3,1x,a,1x,a,2x)',ADVANCE="no",IOSTAT=ios_ith) &
          &'#IC_XYZ  ',PRN_array(i_sat),SVNID,TRIM(BLKTYP),'ICRF'
   WRITE (UNIT=UNIT_IN,FMT='(I5,F19.10)',ADVANCE="no",IOSTAT=ios_ith) INT(orbits_ics_icrf(1,i_sat)), orbits_ics_icrf(2,i_sat)
   WRITE (UNIT=UNIT_IN,FMT= * ,IOSTAT=ios_ith) orbits_ics_icrf(3:Norbits_ics_icrf,i_sat)
! orbitsmatrix_crf(1,3:8,i_sat), ' DR YR BR DC DS YC YS BC BS''(a3,1x,f14.4,f14.6,1x,15(d17.10,1x))'
! IC Kepler
   r_ic = orbits_ics_icrf(3:5,i_sat)
   v_ic = orbits_ics_icrf(6:8,i_sat)
   CALL kepler_z2k(r_ic, v_ic, GM_global, kepler_ic)
   WRITE (UNIT=UNIT_IN,FMT='(a,a3,1x, a,2x)',ADVANCE="no",IOSTAT=ios_ith) &
          &'#IC_Kepler ',PRN_array(i_sat),'Kepler elements &
		  &[Semi-major axis(m) Eccentricity Inclination(deg) Omega_Asc.node(deg) omega_perigee(deg) True-anomaly(deg)]'	  
   WRITE (UNIT=UNIT_IN,FMT= * ,IOSTAT=ios_ith) kepler_ic(1:6)    
END DO 
IF (partials_velocity_cfg > 0) THEN
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO PRN MJD SOD, ICRF [X Y Z ZD YD ZD], ITRF [X Y Z XD YD ZD], &
                                           &Partials [dx/dX dx/dY dx/dZ dx/dXD dx/dYD dx/dZD dY/dX dY/dY dY/dZ &
                                           &dY/dXD dY/dYD dY/dZD ... dx/dRAD1,dx/dRAD1,dx/dRAD3,dx/dRAD3 dx/dRAD4 &
                                           &... dx/dRADN ... dy/dRAD1,dx/dRAD2 ... dy/dRADN ... ]'   
ELSE IF (partials_velocity_cfg == 0) THEN
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO PRN MJD SOD, ICRF [X Y Z ZD YD ZD], ITRF [X Y Z XD YD ZD], &
                                           &Partials [dx/dXo dx/dYo dx/dZo dx/dVxo dx/dVyo dx/dVzo &
                                           &          dy/dXo dy/dYo dy/dZo dy/dVxo dy/dVyo dy/dVzo & 
                                           &          dz/dXo dz/dYo dz/dZo dz/dVxo dz/dVyo dz/dVzo & 
										   &          dx/dRAD1  dx/dRAD2  dx/dRAD3 .....  dx/dRADN &
										   &          dy/dRAD1  dy/dRAD2  dy/dRAD3 .....  dy/dRADN &
										   &          dz/dRAD1  dz/dRAD2  dz/dRAD3 .....  dz/dRADN ]' 
END IF                                                       
WRITE (UNIT=UNIT_IN,FMT= '(A13)' ,IOSTAT=ios_ith) 'End_of_Header'
! ----------------------------------------------------------------------
! Write data to file | Line by line	  
! ----------------------------------------------------------------------
! Write orbit-partials matrix per epoch per satellite
! ----------------------------------------------------------------------
DO i_epoch = 1 , Nepochs
	DO i_sat = 1 , Nsat

!print *,"PRN_array(i_sat)", PRN_array(i_sat)	
!print *,"orbitsmatrix_crf(i_epoch,:,i_sat)", orbitsmatrix_crf(i_epoch,:,i_sat)
	
! Based on the format definition by fmt_wrt 
!WRITE (UNIT=UNIT_IN,FMT=fmt_wrt,IOSTAT=ios_ith) PRN_array(i_sat),' ', orbitsmatrix_crf(i_epoch,:,i_sat)
!WRITE (UNIT=UNIT_IN,FMT= * ,IOSTAT=ios_ith) PRN_array(i_sat),' ', orbitsmatrix_crf(i_epoch,1:8,i_sat), &
!                                            orbitsmatrix_trf(i_epoch,3:8,i_sat),orbitsmatrix_crf(i_epoch,9:NParam,i_sat)
WRITE (UNIT=UNIT_IN,FMT='(A3,2x)',ADVANCE="no",IOSTAT=ios_ith) PRN_array(i_sat)
WRITE (UNIT=UNIT_IN,FMT='(I5)',ADVANCE="no",IOSTAT=ios_ith)  INT(orbitsmatrix_crf(i_epoch,1,i_sat))
WRITE (UNIT=UNIT_IN,FMT='(F19.10)',ADVANCE="no",IOSTAT=ios_ith) orbitsmatrix_crf(i_epoch,2,i_sat)
WRITE (UNIT=UNIT_IN,FMT= * ,IOSTAT=ios_ith) orbitsmatrix_crf(i_epoch,3:8,i_sat), orbitsmatrix_trf(i_epoch,3:8,i_sat), & 
											orbitsmatrix_crf(i_epoch,9:NParam,i_sat)
											
IF (ios_ith /= 0) THEN
   PRINT *, "Error in writing to file: ", TRIM (filename)
   PRINT *, "WRITE IOSTAT=", ios_ith
END IF
	
	END DO
END DO


ENDFILE (UNIT = UNIT_IN) 
CLOSE (UNIT = UNIT_IN)
! ----------------------------------------------------------------------



END SUBROUTINE



End Module

