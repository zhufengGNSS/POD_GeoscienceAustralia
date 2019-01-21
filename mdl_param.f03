MODULE mdl_param

! ---------------------------------------------------------------------------
! Purpose:
!  Module for setting orbit modelling parameters as global variables 
!  to be used from various part of the source code 
! ---------------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	November 2017
! ----------------------------------------------------------------------


      USE mdl_precision
      IMPLICIT NONE
      SAVE 			
	  
	  
! ---------------------------------------------------------------------------
! Time System of orbit ouput files defined by the input configuration file
      CHARACTER (LEN=3) :: TIME_SCALE
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Satellite information file
!      CHARACTER (*), PARAMETER :: satinfoplace ='/data/bitbucket/Tzupang/pod/SATINFO/'
! --------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Orbit arc length (sec)
      REAL (KIND = prec_d) :: orbarc
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Numerical integration method ID
      INTEGER (KIND = prec_int2) :: integmeth	  
! Numerical integration step size (sec)
      REAL (KIND = prec_d) :: integstep
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Initial Conditions
! Initial Epoch: MJD and Seconds since 00h in Terrestrial Time (TT)
      REAL (KIND = prec_d) :: MJD_to
      REAL (KIND = prec_d) :: SEC_to
! Initial State Vector in ICRF
      REAL (KIND = prec_d) :: SVEC_Zo(6)
      REAL (KIND = prec_d) :: SVEC_Zo_ESTIM(6)
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! GNSS satellite PRN number e.g. G03 	  
      CHARACTER (LEN=3) :: PRN
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! External Orbit arrays (prm_orbext.f03, prm_pseudobs.f03)
! Allocatable Arrays	  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbext_kepler, orbext_ICRF, orbext_ITRF
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: pseudobs_ITRF, pseudobs_ICRF
      INTEGER (KIND = prec_int2) :: ORBEXT_glb
! ---------------------------------------------------------------------------


! ---------------------------------------------------------------------------
! Force Model
! ---------------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravitational forces modelling considered (on/off): defined in prm_grav.f03 
INTEGER (KIND = prec_int2) :: FMOD_GRAV(4)
INTEGER (KIND = prec_int2) :: FMOD_GRAVFIELD
INTEGER (KIND = prec_int2) :: FMOD_PLANETS(11)
INTEGER (KIND = prec_int2) :: FMOD_TIDES(5)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity Field Model parameters: setting via module m_gfc.f03
! ----------------------------------------------------------------------
! - GFM_GM:			Earth gravity constant (m^3/sec^2)
! - GFM_ae: 		radius (meters)
! - GFM_degmax:		Maximum degree of the gravity model
! - GFM_tide:		Tide System
! - Cnm, Snm:		Normalized Spherical Harmonic coefficients
! - sCnm, sSnm:		Errors, coefficients variances  
REAL (KIND = prec_q) :: 		GFM_GM
REAL (KIND = prec_q) :: 		GFM_ae
INTEGER (KIND = prec_int8) ::	GFM_degmax 
CHARACTER (LEN=50) ::			GFM_tide 
	  
! Truncation maximum degree and order: user defined (via input parameterization) 
INTEGER (KIND = prec_int8) :: GFM_Nmax, GFM_Mmax

! Spherical Harmonic Coefficients (SHC) lower triangular matrices (Dynamic Allocatable arrays)
REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: GFM_Cnm
REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: GFM_Snm
! Errors/Variances of SHC (Dynamic Allocatable arrays)
REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: GFM_sCnm
REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: GFM_sSnm
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Ocean Tides model
! Truncation maximum degree and order: user defined (via input parameterization) 
INTEGER (KIND = prec_int8) :: OCEAN_Nmax, OCEAN_Mmax
! ---------------------------------------------------------------------------

! ----------------------------------------------------------------------
! Non-gravitational forces modelling considered (on/off): defined in prm_nongrav.f03 
INTEGER (KIND = prec_int2) :: FMOD_NONGRAV(3)
INTEGER (KIND = prec_int2) :: SRP_MOD_glb
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Empirical forces
INTEGER (KIND = prec_int2) :: EMP_param_glb
INTEGER (KIND = prec_int2) :: EMP_Bias_glb(3)
REAL (KIND = prec_q) :: Bias_accel_glb(3), Bias_accel_aposteriori(3)

INTEGER (KIND = prec_int2) :: EMP_CPR_glb(3)
INTEGER (KIND = prec_int2) :: EMP_nCPR_glb
REAL (KIND = prec_q) :: CPR_CS_glb(3,2), CPR_CS_aposteriori(3,2)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Variational Equations
! ----------------------------------------------------------------------
! Number of parameters to be estimated
!INTEGER (KIND = prec_int8) ::	N_PARAM 
INTEGER (KIND = prec_int8) :: NPARAM_glb
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Estimator method
      INTEGER (KIND = prec_int2) :: ESTIM_mode_glb 
      INTEGER (KIND = prec_int2) :: ESTIM_iter_glb 
! ----------------------------------------------------------------------

END