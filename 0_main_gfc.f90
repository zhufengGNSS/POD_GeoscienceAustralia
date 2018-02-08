PROGRAM main_gfc


! ----------------------------------------------------------------------
! Program: main_gfc.f90
! ----------------------------------------------------------------------
! Purpose:
!  Read global gravity field model data in ICGEM format
! ----------------------------------------------------------------------
! Called subroutines:
!  gfc1.f90
!  gfc2.f90
!  gfc3_iers.f90
!  iau_CAL2JD.for 
!  write_array.f90
! ----------------------------------------------------------------------
! Used modules: 
!  mdl_precision.f90
!  mdl_num.f90
!  mdl_gfc.f90
!  mdl_write.f90
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia           September 2015
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_gfc
      USE mdl_write
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: n_max, m_max, Ntrunc 
      INTEGER (KIND = prec_int2) :: sigma_shc
      CHARACTER (LEN=300) :: gfmfilename, filename	    
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! INPUT 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity Field model
! ----------------------------------------------------------------------
      !gfmfilename = 'egm2008.gfc'
      !gfmfilename = "go_cons_gcf_2_dir_r4.gfc"
      gfmfilename = 'goco05s.gfc'
! ----------------------------------------------------------------------
! Degree truncation limit (Spherical harmonics expansion series)
	  Ntrunc = 5
! Errors of spherical harmonic coefficients: Store errors in arrays or Not
      sigma_shc = 0  
! ----------------------------------------------------------------------
      PRINT *, "gfmfilename: ", TRIM(gfmfilename)
      print *,"Truncation Nmax:", Ntrunc 
 
! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Gravity Field model: Spherical Harmonic Coefficients (icgem format) 
! ----------------------------------------------------------------------
      CALL gfc1 (gfmfilename,Ntrunc,sigma_shc)
! ----------------------------------------------------------------------
      print *,"Cnm=", Cnm
      print *,"Snm=", Snm

	  
! ----------------------------------------------------------------------
! GFM parameters : mdl_gfc.f90
! ----------------------------------------------------------------------
      GM_global = GM_gfc
      IF (Ntrunc == -1) THEN
         n_max = Nmax_gfc
      ELSE
         n_max = Ntrunc
      END IF  
      m_max = n_max
      print *,"Nmax_gfc=", Nmax_gfc, n_max, m_max 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Write Cnm,Snm lower triangulation matrices to ascii files
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
