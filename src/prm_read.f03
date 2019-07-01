SUBROUTINE prm_read (PRMfname)


! ----------------------------------------------------------------------
! SUBROUTINE: prm_read.f03
! ----------------------------------------------------------------------
! Purpose:
!  Read the orbit parameterization based on the input configuration file 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Cooperative Research Centre for Spatial Information, Australia
! Created:	04 October 2017
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_eop
      USE mdl_arr
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN) :: PRMfname				

! OUT

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      CHARACTER (LEN=50) :: REF_Zo, time_system, Sat_PRN
      INTEGER (KIND = prec_int4) :: time_in
      CHARACTER (LEN=1) :: GNSSid
      INTEGER (KIND = 4) :: PRN_sp3
      INTEGER (KIND = 4) :: PRN_eclips

      CHARACTER (LEN=7) :: integrator

      REAL (KIND = prec_d) :: Zo(6), ro(3), vo(3)	  
      INTEGER IY, IM, ID, J_flag
      DOUBLE PRECISION DJM0, Sec, FD
      REAL (KIND = prec_d) :: mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC

      REAL (KIND = prec_d) :: r_TRS(3), v_TRS(3)
      REAL (KIND = prec_d) :: r_CRS(3), v_CRS(3), v_CRS_1(3), v_CRS_2(3)

      REAL (KIND = prec_d) :: EOP_cr(7)
      REAL (KIND = prec_d) :: CRS2TRS(3,3), TRS2CRS(3,3)
      REAL (KIND = prec_d) :: d_CRS2TRS(3,3), d_TRS2CRS(3,3)

      REAL (KIND = prec_d) :: GM

      CHARACTER (LEN=100) :: fmt_line
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, read_i
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios
      INTEGER (KIND = prec_int2) :: ios_line, ios_key, ios_data
      INTEGER (KIND = prec_int2) :: space_i
      INTEGER (KIND = prec_int2) :: AllocateStatus
      CHARACTER (LEN=7) :: Format1, Format2, Format3
      CHARACTER (LEN=500) :: line_ith	  
      CHARACTER (LEN=150) :: word1_ln, word_i, t0	  
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: ForceMod
! Gravity Field model variables
      INTEGER (KIND = prec_int8) :: Ntrunc
      INTEGER (KIND = prec_int8) :: Nmax 
      INTEGER (KIND = prec_int2) :: gfc_opt
      CHARACTER (LEN=300) :: gfmfilename
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Orbit parameterization INPUT file read:
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      UNIT_IN = 9  												
      Format1 = '(A)'
      Format2 = '(F)'
      Format3 = '(I100)'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Open .in file
      OPEN (UNIT = UNIT_IN, FILE = TRIM (PRMfname), IOSTAT = ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", PRMfname
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Read input file
i = 0
DO

READ (UNIT=UNIT_IN,FMT=Format1,IOSTAT=ios_line) line_ith
i = i + 1
! PRINT *, "READ Line (i,ios):", i, ios_line

! ----------------------------------------------------------------------
! End of file
         IF (ios_line < 0) THEN
!            PRINT *, "End of file, i=", i
            EXIT		
         END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1st Word of Line ith
READ (line_ith, * , IOSTAT=ios_data) word1_ln  ! 1st word
!READ (line_ith, * , IOSTAT=ios_data) word1_ln, charN 
! ----------------------------------------------------------------------
!PRINT *, "word1_ln: ", word1_ln


! ----------------------------------------------------------------------
! Parameters Keywords read 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Reference Frame of the input initial state vector
! 1. ITRF
! 2. ICRF
! 3. Kepler Elements (in ICRF)
IF (word1_ln == "Reference_frame") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, REF_Zo 
END IF
! ----------------------------------------------------------------------
		 
! ----------------------------------------------------------------------
! Time System of the input initial epoch:
! 1. TT
! 2. GPS time
! 3. UTC
! 4. TAI
IF (word1_ln == "Time_scale") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, time_system 
END IF

If (time_system == 'TT') Then
	time_in = 1
Else if (time_system == 'GPS') then
	time_in = 2
Else if (time_system == 'UTC') then
	time_in = 3
Else if (time_system == 'TAI') then
	time_in = 4
End If	

! TIME_System: Global variable in module mdl_param.f03
! defining the time scale of the ouput files of the computed orbits
TIME_SCALE  = time_system	 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! GNSS Satellite PRN number 
! ----------------------------------------------------------------------
! GNSS PRN number: Constellation ID + Number e.g. G04
! ----------------------------------------------------------------------
! GNSSid:	GNSS constellation id letter
! 			G: GPS
! 			R: GLONASS
! 			E: Galileo
! 			C: BDS (BeiDou)
! 			J: QZSS
! 			S: SBAS
!
! PRN_sp3:	PRN numbering as adopted in the sp3 format (Numerical value following the GNSS constellation id letter)
! ----------------------------------------------------------------------
IF (word1_ln == "Satellite_PRN") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Sat_PRN 
END IF
PRN = TRIM(Sat_PRN)

!fmt_line = '(A1,I2.2)'
!READ (PRN, fmt_line , IOSTAT=ios) GNSSid, PRN_sp3
!print *, "GNSSid ", GNSSid
!print *, "PRN_sp3 ", PRN_sp3
! ----------------------------------------------------------------------
		 
! ----------------------------------------------------------------------
! Orbit arc length (sec)
IF (word1_ln == "Orbit_arc_length") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, orbarc 
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Initial Conditions

! Intial Epoch
IF (word1_ln == "Year") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, IY 
END IF

IF (word1_ln == "Month") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, IM 
END IF

IF (word1_ln == "Day") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ID 
END IF

IF (word1_ln == "Seconds") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Sec 
END IF

! Initial State vector (m)      
IF (word1_ln == "state_vector") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Zo 
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Earth Orientation                                                       
! ----------------------------------------------------------------------
! EOP settings are saved as global variables in the module mdl_EOP.f90 (to be used by the subroutine EOP.f90)

! EOP data solution options
! 1. IERS C04
! 2. IERS RS/PC Daily (finals2000A.daily)
! 3. IGS ultra-rapid ERP + IERS RS/PC Daily (dX,dY)
IF (word1_ln == "EOP_data_sol") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, EOP_sol 
END IF

! EOP filename
IF (word1_ln == "EOP_filename") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, EOP_fname 
END IF

! ERP filename
IF (word1_ln == "ERP_filename") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ERP_fname 
END IF

! EOP interpolation points
IF (word1_ln == "EOP_interpolation_points") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, EOP_Nint 
END IF

! IAU Precession-Nutation model:
! 2. IAU2000A:			iau_model = 2000
! 3. IAU2006/2000A:		iau_model = 2006
IF (word1_ln == "iau_pn_model") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, iau_model 
END IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Numerical Integrator                                                    
! ----------------------------------------------------------------------
! Numerical Integration methods
! 1. RKN7(6)8:	Runge-Kutta-Nystrom 7th order   
! 2. RK4:		Runge-Kutta 4th order
! 3. RK8(7)13:	Runge-Kutta 8th order
IF (word1_ln == "integrator_meth") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, integrator 
END IF

If (integrator =='RKN7') Then
	integmeth = 1
Else if (integrator =='RK4') then
	integmeth = 2
Else if (integrator =='RK8') then
	integmeth = 3
End IF

! Numerical Integration Stepsize (in sec)
IF (word1_ln == "integrator_step") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, integstep 
END IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Gravitational forces considered: 
! ----------------------------------------------------------------------
! FMOD_GRAV setting (global variable in the module mdl_param.f03)
! 
! FMOD_GRAV(1) : Earth gravity Field
! FMOD_GRAV(2) : Planetary/Lunar orbital perturbations
! FMOD_GRAV(3) : Tides
! FMOD_GRAV(4) : Relativistic Effects
! 
! 1. FMOD_GRAV(i)=1 : Effect is considered
! 2. FMOD_GRAV(i)=0 : Effect is neglected 
!
! ----------------------------------------------------------------------
! Gravity Field
IF (word1_ln == "Gravity_field") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_GRAV(1) = ForceMod
	!print *,"ForceMod", ForceMod
END IF

! Planetary/Lunar Perturbations
IF (word1_ln == "Planets_perturbations") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_GRAV(2) = ForceMod
	!print *,"ForceMod", ForceMod
END IF

! Tidal effects
IF (word1_ln == "Tidal_effects") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_GRAV(3) = ForceMod
	!print *,"ForceMod", ForceMod
END IF

! Relativistic effects
IF (word1_ln == "Relativistic_effects") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_GRAV(4) = ForceMod
	!print *,"ForceMod", ForceMod
END IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Earth Gravity Field
! ----------------------------------------------------------------------
! Gravity Field model options:
! 0. Central gravity field (Keplerian orbit only)	: gravity_model = 0
! 1. Static global gravity field model     		    : gravity_model = 1
! 2. Time-variable global gravity field model 		: gravity_model = 2
! 3. IERS conventional geopotential model 			: gravity_model = 3
IF (word1_ln == "gravity_model") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, gfc_opt 
END IF

! FMOD_GRAVFIELD global variable in mdl_param.f03
IF (gfc_opt == 0) Then
	FMOD_GRAVFIELD = 0
ELSE
	FMOD_GRAVFIELD = 1
END IF

! Gravity Field model
IF (word1_ln == "gravity_model_filename") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, gfmfilename 
END IF

! Degree truncation limit (Spherical harmonics expansion series)
IF (word1_ln == "degree_max") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Nmax 
	Ntrunc = Nmax
END IF

! Global variables set in module mdl_param.f03
GFM_Nmax = Nmax
GFM_Mmax = Nmax

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------	
! Tides
! ----------------------------------------------------------------------	
! FMOD_TIDES setting (global variable in the module mdl_param.f03)
! 
! FMOD_TIDES(1) : Solid Earth Tides: Frequency indepedent terms
! FMOD_TIDES(2) : Solid Earth Tides: Frequency depedent terms
! FMOD_TIDES(3) : Ocean Tides
! FMOD_TIDES(4) : Solid Earth Pole Tide 
! FMOD_TIDES(5) : Ocean Pole Tide Tide
! 
! 1. FMOD_TIDES(i)=1 : Effect is considered
! 2. FMOD_TIDES(i)=0 : Effect is neglected 
!
! FMOD_TIDES = (/ 1, 1, 1, 1, 1 /)

IF (word1_ln == "solid_tides_nonfreq") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_TIDES(1) = ForceMod
END IF

IF (word1_ln == "solid_tides_freq") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_TIDES(2) = ForceMod
END IF

IF (word1_ln == "ocean_tides") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_TIDES(3) = ForceMod
END IF

IF (word1_ln == "solid_earth_pole_tide") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_TIDES(4) = ForceMod
END IF

IF (word1_ln == "ocean_pole_tide") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_TIDES(5) = ForceMod
END IF
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Non-gravitational forces considered: 
! ----------------------------------------------------------------------
! FMOD_NONGRAV setting (global variable in the module mdl_param.f03)
! 
! FMOD_NONGRAV(1) : Solar Radiation
! FMOD_NONGRAV(2) : Earth Radiation
! FMOD_NONGRAV(3) : Antenna Thrust
! 
! 1. FMOD_NONGRAV(i)=1 : Effect is considered
! 2. FMOD_NONGRAV(i)=0 : Effect is neglected 
!
! ----------------------------------------------------------------------
! Solar Radiation
IF (word1_ln == "Solar_radiation") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_NONGRAV(1) = ForceMod
	!print *,"ForceMod", ForceMod
END IF

! Earth Radiation
IF (word1_ln == "Earth_radiation") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_NONGRAV(2) = ForceMod
	!print *,"ForceMod", ForceMod
END IF

! Antenna Thrust
IF (word1_ln == "Antenna_thrust") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ForceMod 
	FMOD_NONGRAV(3) = ForceMod
	!print *,"ForceMod", ForceMod
END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Solar Radiation
! ----------------------------------------------------------------------
! Solar Radiation pressure model options:
! 1. Cannonball model
! 2. Box-wing model
! 3. ECOM model
IF (word1_ln == "SRP_model") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, SRP_MOD_glb 
END IF
! ----------------------------------------------------------------------



END DO
CLOSE (UNIT=UNIT_IN)
! Close of input parameterization file
! ----------------------------------------------------------------------





! ----------------------------------------------------------------------
! Initial Conditions: Time & Reference frame transformation
! ----------------------------------------------------------------------
! MJD of input initial epoch including fraction of the day	  
CALL iau_CAL2JD ( IY, IM, ID, DJM0, mjd, J_flag )
mjd = mjd + Sec / 86400.D0	  

If (time_in == 1) Then

! ----------------------------------------------------------------------
! Initial Epoch in TT
MJD_to = mjd
! Seconds since 00h
SEC_to = Sec
! ----------------------------------------------------------------------

Else

! ----------------------------------------------------------------------
! "Time Systems" transformation											 
      IF (time_in == 1) THEN
	     CALL time_TT (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)
      ELSE IF (time_in == 2) THEN 
	     CALL time_GPS (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)
      ELSE IF (time_in == 3) THEN 
	     CALL time_UTC (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)
      ELSE IF (time_in == 4) THEN 
         CALL time_TAI (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)
      END IF 
! ----------------------------------------------------------------------
! GPS week and Seconds of GPS week
! CALL time_GPSweek (mjd_GPS , GPS_week, GPS_wsec, GPSweek_mod1024)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Initial Epoch in TT
MJD_to = mjd_TT
! Seconds since 00h
SEC_to = (MJD_to - INT(MJD_to)) * (24.D0 * 3600.D0)
! ----------------------------------------------------------------------

End If


! ----------------------------------------------------------------------
! Reference frame of initial state vector
! ----------------------------------------------------------------------
If (REF_Zo == 'ICRF') Then

    SVEC_Zo = Zo
	
Else If (REF_Zo == 'ITRF') Then

r_TRS(1:3) = Zo(1:3)  
v_TRS(1:3) = Zo(4:6)  

! ICRF-ITRF transformation matrix (including derivatives) based on EOP data
CALL EOP (mjd_TT, EOP_cr, CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS)

! ITRF to ICRF 
! r_CRS = TRS2CRS * r_TRS
      CALL matrix_Rr (TRS2CRS,r_TRS , r_CRS)	  
! v_CRS = TRS2CRS * v_TRS + d_TRS2CRS * r_TRS
      CALL matrix_Rr (TRS2CRS,v_TRS , v_CRS_1)
      CALL matrix_Rr (d_TRS2CRS,r_TRS , v_CRS_2)
      v_CRS = v_CRS_1 + v_CRS_2
	  
	  SVEC_Zo(1:3) = r_CRS(1:3)	
	  SVEC_Zo(4:6) = v_CRS(1:3)	

Else If (REF_Zo == 'Kepler') Then

! GM global value via module mdl_num.f90
	 GM = GM_global

	 !Call kepler_k2z (kepler, GM , r,v)
	 Call kepler_k2z (Zo, GM , ro, vo)
	 SVEC_Zo(1:3) = ro
	 SVEC_Zo(4:6) = vo

End If
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Forces model settings
! ----------------------------------------------------------------------
! Gravitational Effects
!Call prm_grav (PRMfname)

! Non-gravitational Effects
!Call prm_nongrav (PRMfname)

! Empirical parameters/accelerations

! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Observation model                                                       
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------







if (1<0) then
! ----------------------------------------------------------------------
PRINT *,"--------------------- INPUT ----------------------------"
PRINT *, "Reference Frame:               ", REF_Zo
PRINT *, "Time System:                   ", time_system
PRINT *, "GNSS Satellite PRN:            ", PRN
PRINT *, "Orbit arc length:              ", orbarc
PRINT *, "Epoch:                         ", IY, IM, ID, Sec
PRINT *, "ro:                 		     ", Zo(1:3)
PRINT *, "vo:                 		     ", Zo(4:6)

PRINT *, "EOP data solution ID:          ", EOP_sol
PRINT *, "EOP data file name:            ", TRIM(EOP_fname)
PRINT *, "ERP data file name:            ", TRIM(ERP_fname) 
PRINT *, "EOP interpolation points       ", EOP_Nint
PRINT *, "IAU Precession-Nutation model: ", iau_model

PRINT *, "Numerical Integrator method    ", integrator
PRINT *, "Numerical Integrator ID        ", integmeth
PRINT *, "Numerical Integration step     ", integstep
PRINT *,"--------------------------------------------------------"
PRINT *," "
! ----------------------------------------------------------------------
end if




END
