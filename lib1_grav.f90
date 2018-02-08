PROGRAM lib1_grav


! ----------------------------------------------------------------------
! Program: lib1_grav.f90
! ----------------------------------------------------------------------
! Purpose:
!  Gravitational components due to the Earth gravity field
! ----------------------------------------------------------------------
! Called subroutines:
!  gfc1.f90
!  gfc2.f90
!  gfc3_iers.f90
!  force_gfm.f90
!  iau_CAL2JD.for 
!  write_array.f90
! ----------------------------------------------------------------------
! Used modules: 
!  mdl_precision.f90
!  mdl_num.f90
!  mdl_legendre.f90
!  mdl_legendre1.f90
!  mdl_gfc.f90
!  mdl_write.f90
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia           September 2015
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_legendre
      USE mdl_legendre1
      USE mdl_gfc
      USE mdl_write
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: x,y,z
      REAL (KIND = prec_q), DIMENSION(3) :: r 
! ----------------------------------------------------------------------
! Gravity Field model variables
      INTEGER (KIND = prec_int8) :: n_max, m_max, Ntrunc, Ntrunc_tv 
      INTEGER (KIND = prec_int2) :: sigma_shc, gfc_opt
      CHARACTER (LEN=300) :: gfmfilename, filename	    
! ----------------------------------------------------------------------
      INTEGER IY, IM, ID, J_flag
      DOUBLE PRECISION DJM0, sec, mjd_t
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: fx,fy,fz, fg
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! INPUT 
! ----------------------------------------------------------------------
! Position vector (GPS)  
      x = 17835.809256D0 * 1.0D3
      y =  4999.483252D0 * 1.0D3
      z = 19053.904501D0 * 1.0D3

	  r(1) = x
	  r(2) = y
	  r(3) = z
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Earth Gravity Field
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1. Global gravity field model (static model)    		: set gfc_opt = 1
! 2. Global gravity field model (time-variable model)	: set gfc_opt = 2
! 2. IERS conventional geopotential model 				: set gfc_opt = 3
      gfc_opt = 2
      PRINT *,"gfc_opt", gfc_opt
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity Field model
      !gfmfilename = 'egm2008.gfc'   ! Static Geopotential adopted by IERS Conv. 2010
      !gfmfilename = "go_cons_gcf_2_dir_r4.gfc"
      gfmfilename = 'goco05s.gfc'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Degree truncation limit (Spherical harmonics expansion series)
	  Ntrunc = 20
! Degree truncation limit for the time-variable coefficients (Only when a time-variable gravity model is selected and gfc_opt=2)
      Ntrunc_tv = 0
! Errors of spherical harmonic coefficients: Store errors in arrays or Not
      sigma_shc = 0  
! ----------------------------------------------------------------------
      PRINT *, "gfmfilename: ", TRIM(gfmfilename)
      print *,"Truncation Nmax, Nmax_tv:", Ntrunc, Ntrunc_tv 


	  
! ----------------------------------------------------------------------
! Epoch | Required for the time-variable gravity field effects
      IY = 2015
      IM = 7
      ID = 22
      sec = 600.0D0  
      CALL iau_CAL2JD ( IY, IM, ID, DJM0, mjd_t, J_flag )
      mjd_t = mjd_t + sec / 86400D0
! ----------------------------------------------------------------------

  
! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Gravity Field model: Spherical Harmonic Coefficients (icgem format) 
! ----------------------------------------------------------------------
      if (gfc_opt == 1) then 
          CALL gfc1 (gfmfilename,Ntrunc,sigma_shc)
! ----------------------------------------------------------------------
      else if (gfc_opt == 2) then 	  
          CALL gfc2 (gfmfilename,Ntrunc,sigma_shc, mjd_t, Ntrunc_tv) 
! ----------------------------------------------------------------------
      else if (gfc_opt == 2) then 	  
          CALL gfc3_iers (gfmfilename, Ntrunc, sigma_shc, mjd_t)
! ----------------------------------------------------------------------
      end if
      !print *,"Cnm=", Cnm
      !print *,"Snm=", Snm

	  
! ----------------------------------------------------------------------
! GFM parameters : mdl_gfc.f90
! ----------------------------------------------------------------------
      !GM_global = GM_gfc
      IF (Ntrunc == -1) THEN
         n_max = Nmax_gfc
      ELSE
         n_max = Ntrunc
      END IF  
      m_max = n_max
! ----------------------------------------------------------------------
      print *,"Nmax_gfc=", Nmax_gfc, n_max, m_max 



! ----------------------------------------------------------------------
! Gravitational cartesian components
      CALL force_gfm (r,n_max,m_max , fx,fy,fz) 		
      fg = sqrt(fx**2 + fy**2 + fz**2)
! ----------------------------------------------------------------------
      print *,"fx=", fx
      print *,"fy=", fy
      print *,"fz=", fz
      print *,"fg=", fg	  



! ----------------------------------------------------------------------
! Write Cnm and Snm matrices to ascii files
! ----------------------------------------------------------------------
! Cnm
      ALLOCATE (wrtArray(n_max+1,n_max+1), STAT=AllocateStatus)
      wrtArray = Cnm
      filename = "Cnm.wrt"
      CALL write_array (filename)
      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
! Snm
      ALLOCATE (wrtArray(n_max+1,n_max+1), STAT=AllocateStatus)
      wrtArray = Snm
      filename = "Snm.wrt"
      CALL write_array (filename)
      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Arrays DeALLOCATION
! ----------------------------------------------------------------------
      DEALLOCATE (Cnm, STAT = DeAllocateStatus)
!      print *,"Cnm DeAllocateStatus=", DeAllocateStatus
      DEALLOCATE (Snm, STAT = DeAllocateStatus)
!      print *,"Snm DeAllocateStatus=", DeAllocateStatus
! ----------------------------------------------------------------------


END
