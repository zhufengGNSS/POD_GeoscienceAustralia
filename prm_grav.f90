SUBROUTINE prm_grav ( )


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
      USE mdl_legendre
      USE mdl_legendre1
      USE mdl_gfc
      USE mdl_planets
      USE mdl_tides	  
      USE mdl_write
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
! Gravity Field model variables
      INTEGER (KIND = prec_int8) :: Ntrunc, Ntrunc_tv 
      !INTEGER (KIND = prec_int8) :: n_max, m_max
      INTEGER (KIND = prec_int2) :: sigma_shc, gfc_opt
      CHARACTER (LEN=300) :: gfmfilename	   
! ----------------------------------------------------------------------
! Planetary/Lunar orbits variables
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out
! ----------------------------------------------------------------------
! Tides
      CHARACTER (LEN=100) :: FESxxfname				
      INTEGER (KIND = prec_int8) :: ocean_Nmax, ocean_Mmax
      INTEGER (KIND = prec_int2) :: tide_system
! ----------------------------------------------------------------------
      CHARACTER (LEN=300) :: filename	   
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Earth Gravity Field
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity field model type options
! 1. Global gravity field model (static model)    		: set gfc_opt = 1
! 2. Global gravity field model (time-variable model)	: set gfc_opt = 2
! 2. IERS conventional geopotential model 				: set gfc_opt = 3
      gfc_opt = 2
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity Field model
      gfmfilename = 'goco05s.gfc'
      !gfmfilename = "go_cons_gcf_2_dir_r4.gfc"
      !gfmfilename = 'egm2008.gfc'   
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Degree truncation limit (Spherical harmonics expansion series)
	  Ntrunc = 20
! Degree truncation limit for the time-variable coefficients (Only when a time-variable gravity model is selected and gfc_opt=2)
      Ntrunc_tv = 0
! Errors of spherical harmonic coefficients: Store errors in arrays or Not
      sigma_shc = 0  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity model's data: Spherical Harmonic Coefficients (icgem format)
! ----------------------------------------------------------------------
! Read gravity model parameters and data
! Store the model parameters as global variables via module mdl_gfc.f90 
! Store the spherical harmonic coefficients in dynamic allocatable arrays via module mdl_gfc.f90 
      if (gfc_opt == 1) then 
          CALL gfc1 (gfmfilename,Ntrunc,sigma_shc)
      else if (gfc_opt == 2) then 	  
          CALL gfc2 (gfmfilename,Ntrunc,sigma_shc, mjd_t, Ntrunc_tv) 
      else if (gfc_opt == 2) then 	  
          CALL gfc3_iers (gfmfilename, Ntrunc, sigma_shc, mjd_t)
      end if
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Set maximum degree and order 	  
	! Nmax
      IF (Ntrunc == -1) THEN
         Nmax = Nmax_gfc
      ELSE
         Nmax = Ntrunc
      END IF  
	! Mmax
      Mmax = Nmax
! ----------------------------------------------------------------------

! End of parameters setting for Gravity Field
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Planetary/Lunar orbit data
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Planetary/Lunar precise ephemeris DE
      fname_out = 'DE.430' 
      fname_header = 'header.430_229'
      fname_data = 'ascp1950.430'
      CALL CATfile (fname_header,fname_data,fname_out)

! DE ephemeris data processing
! Store selected DE data to global variables via the module mdl_planets.f90 
      CALL asc2eph (fname_out)

! GM gravity constants of the solar system bodies ! Input/Output arguments via the module mdl_planets.f90
      CALL GM_de
! ----------------------------------------------------------------------

! End of parameters setting for Planetary/Lunar orbit data
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
!  Tidal effects to orbits: Solid Earth tides, Ocean Tides and Pole Tide
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Tide System (Tide-free or Zero-tide)
! ----------------------------------------------------------------------
! Tide system is defined by the input gravity field model (tide_gfc variable)
! tide_gfc variable is set in the module mdl_gfc.f90
! 1. Tide-free: tide_system = 0
! 2. Zero-tide: tide_system = 1
if (tide_gfc == 'tide_free') then
      tide_system = 0 
else if (tide_gfc == 'zero_tide') then
      tide_system = 1 
end if 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Ocean Tides model
! ----------------------------------------------------------------------
	  FESxxfname = 'fes2004_Cnm-Snm.dat'

! Read ocean tides model data: Spherical harmonic coefficients corrections
! The spherical harmoncis corrections are stored in dynamic allocatable arrays through the module mdl_tides.f90
      CALL tides_fes2004(FESxxfname)
	  
! Maximum degree and order of ocean tides model spherical harmonic coefficients corrections	  
      ocean_Nmax = Nmax
	  ocean_Mmax = ocean_Nmax
! ----------------------------------------------------------------------


! End of parameters setting for Tidal effects
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Relativity parameters: beta and gama
! ----------------------------------------------------------------------
! PPN (parameterized post-Newtonian) parameters  
! PPN parameters are set equal to 1 in General Relativity

! beta and gama parameters are declared as global variables in the module mdl_num.f90
! ----------------------------------------------------------------------




if (1<0) then
      PRINT *,"gfc_opt", gfc_opt
      PRINT *, "gfmfilename: ", TRIM(gfmfilename)
      print *,"Truncation Nmax, Nmax_tv:", Ntrunc, Ntrunc_tv 
      !print *,"Cnm=", Cnm
      !print *,"Snm=", Snm
      print *,"Nmax_gfc=", Nmax_gfc, Nmax, m_max 
end if



if (1<0) then
! ----------------------------------------------------------------------
! Write spherical harmonic coefficients, Cnm and Snm matrices, to ascii files
! ----------------------------------------------------------------------
! Cnm
      ALLOCATE (wrtArray(Nmax+1,Nmax+1), STAT=AllocateStatus)
      wrtArray = Cnm
      filename = "Cnm.out"
      CALL write_array (filename)
      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
! Snm
      ALLOCATE (wrtArray(Nmax+1,Nmax+1), STAT=AllocateStatus)
      wrtArray = Snm
      filename = "Snm.out"
      CALL write_array (filename)
      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
end if




END
