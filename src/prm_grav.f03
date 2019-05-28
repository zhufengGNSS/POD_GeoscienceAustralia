SUBROUTINE prm_grav (PRMfname)


! ----------------------------------------------------------------------
! SUBROUTINE: prm_grav.f90
! ----------------------------------------------------------------------
! Purpose:
!  Read/Set the parameterization regarding the gravitational effects to orbit
!  Read the data and set the parameterization to global variables through modules 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Cooperative Research Centre for Spatial Information, Australia
! Created:	20 September 2017
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      !USE mdl_gfc
      !USE m_gfc
      !USE m_gfc3
      !USE mdl_legendre
      !USE mdl_legendre1
      !USE mdl_planets
      !USE mdl_tides	  
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN) :: PRMfname				
! OUT

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: ForceMod
	  
! Gravity Field model variables
      INTEGER (KIND = prec_int8) :: Ntrunc, Ntrunc_tv 
      INTEGER (KIND = prec_int8) :: Nmax, Mmax 
      !INTEGER (KIND = prec_int8) :: n_max, m_max
      INTEGER (KIND = prec_int2) :: sigma_shc, gfc_opt
      CHARACTER (LEN=300) :: gfmfilename
      REAL (KIND = prec_q) :: mjd_t

      REAL (KIND = prec_q) :: GM_gfc, ae_gfc
      INTEGER (KIND = prec_int8) :: Nmax_gfc 
      CHARACTER (LEN=50) :: tide_gfc  
	  
!! Spherical Harmonic Coefficients (SHC) lower triangular matrices (Dynamic Allocatable arrays)
!      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: Cnm
!      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: Snm
! Errors/Variances of SHC (Dynamic Allocatable arrays)
!      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: sCnm
!      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: sSnm  
! ----------------------------------------------------------------------
! Planetary/Lunar orbits variables
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out
! ----------------------------------------------------------------------
! Tides
      CHARACTER (LEN=100) :: FESxxfname				
!      INTEGER (KIND = prec_int8) :: ocean_Nmax, ocean_Mmax
!      INTEGER (KIND = prec_int2) :: tide_system
! ----------------------------------------------------------------------
      CHARACTER (LEN=300) :: filename	   
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, read_i
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios
      INTEGER (KIND = prec_int2) :: ios_line, ios_key, ios_data
      INTEGER (KIND = prec_int2) :: space_i
      CHARACTER (LEN=7) :: Format1, Format2, Format3
      CHARACTER (LEN=500) :: line_ith	  
      CHARACTER (LEN=150) :: word1_ln, word_i, t0	  
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
!PRINT *, "READ Line (i,ios):", i, ios_line

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
! gfc_opt = 2
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
! gfmfilename = 'goco05s.gfc'
END IF

! Degree truncation limit (Spherical harmonics expansion series)
IF (word1_ln == "degree_max") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Nmax 
	Ntrunc = Nmax
END IF

! Degree truncation limit for the time-variable coefficients (case: gfc_opt=2)
IF (word1_ln == "degree_max_timevar") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Nmax 
	Ntrunc_tv = Nmax
END IF
	  
! Errors of spherical harmonic coefficients: Store errors in arrays or Not
sigma_shc = 0  
! ----------------------------------------------------------------------	


! ----------------------------------------------------------------------	
! Planetary/Lunar precise ephemeris DE
! ----------------------------------------------------------------------	
! DE Header data file
IF (word1_ln == "DE_fname_header") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, fname_header 
!fname_header = 'header.430_229'
END IF

! DE Header data file
IF (word1_ln == "DE_fname_data") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, fname_data 
!fname_data = 'ascp1950.430'
END IF
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
! Ocean Tides model
! ----------------------------------------------------------------------
! Model file name
IF (word1_ln == "ocean_tides_model_file") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, filename 
	FESxxfname = TRIM(filename)
END IF

! Maximum degree and order of ocean tides model (spherical harmonic coefficients corrections)	  
IF (word1_ln == "ocean_tides_model_deg") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Nmax 
	OCEAN_Nmax = Nmax
	OCEAN_Mmax = Nmax
END IF
! ----------------------------------------------------------------------


END DO
CLOSE (UNIT=UNIT_IN)
! Close of input parameterization file
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Relativity parameters: beta and gama
! ----------------------------------------------------------------------
! PPN (parameterized post-Newtonian) parameters  
! PPN parameters are set equal to 1 in General Relativity

! beta and gama parameters are declared as global variables in the module mdl_num.f90
! ----------------------------------------------------------------------



END
