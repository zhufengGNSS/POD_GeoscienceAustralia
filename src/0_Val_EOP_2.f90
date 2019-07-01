      program Val_EOP_2


! ----------------------------------------------------------------------
! Program:	Val_EOP_2.f90
! ----------------------------------------------------------------------
! Purpose:
!  Validation test of TRS - EOP data and Earth Orientation Matrix:
!  Extended numerical comparison between TRS (lib7_trs.f90) and IERS web calculation service
!  in terms of the Earth Orientation Matrix (EOM) elements
! ----------------------------------------------------------------------
! Called subroutines (basic):
! - time_TT.f90, time_TAI.f90, time_GPS.f90, time_UTC.f90 : Time scales change
! - eop_cor.f90:	IERS EOP data processing and tidal corrections (eop_interp.f90, interp.f)
! - eop_urapid.f90 (eop_igu_int.f90) : IGS ultra-rapid ERP data processing
! - CRS_TRS.f90:	ICRF-ITRF transformation matrix (direct/inverse & derivatives)
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia               April 2016
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_write
      USE mdl_arr
      IMPLICIT NONE

! ----------------------------------------------------------------------
      DOUBLE PRECISION arcsec2rad
      INTEGER (KIND = prec_int4) :: time_in
      INTEGER IY, IM, ID, J_flag
      DOUBLE PRECISION DJM0, sec, FD
      REAL (KIND = prec_d) :: mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC, mjd_UT1
      CHARACTER (LEN=50) :: EOP_fname
      INTEGER (KIND = prec_int4) :: n_interp
      INTEGER (KIND = prec_int8) :: mjd_UTC_day
      INTEGER (KIND = prec_int1) :: EOP_sol
      INTEGER (KIND = prec_int2) :: iau_model
      DOUBLE PRECISION EOP_data(7), EOP_cr(7), delta_EOP(10), EOP_cr2(7)
      DOUBLE PRECISION GCRS2ITRS(3,3), ITRS2GCRS(3,3)
      DOUBLE PRECISION GCRS2ITRS_2(3,3), ITRS2GCRS_2(3,3), dTM(3,3)
      DOUBLE PRECISION CRS2TRS(3,3), TRS2CRS(3,3), d_CRS2TRS(3,3), d_TRS2CRS(3,3)
! ----------------------------------------------------------------------
      !REAL (KIND = prec_d) :: dif_MJD, mjd_iers, TAIsec_iers, TAI_UTC_iers
      !DOUBLE PRECISION dif_MJD, mjd_iers, TAIsec_iers, TAI_UTC_iers
      !DOUBLE PRECISION dif_EOM(3,3), R_iers(3,3)
      !DOUBLE PRECISION dif_PNM(3,3), PN_iers(3,3)
      !DOUBLE PRECISION dif_X, dif_Y, X_iersweb, Y_iersweb
! ----------------------------------------------------------------------
      CHARACTER (LEN=50) :: ERP_fname
      REAL (KIND = prec_d) :: mjd_igu
! ----------------------------------------------------------------------
! Allocatable Arrays	  
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: R_lib, R_iers, d_eom 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int4) :: epoch_i, t_step, t_period, Nepochs, sz1, sz2, i,j,i2, iers_web
	  CHARACTER (LEN=50) :: fname_eom_iers, fname_eom_iers_ar, fname_eom_lib, fname_dEOM 
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
      REAL (KIND = prec_q) :: stat_q(5)




! ----------------------------------------------------------------------
! PI precision quad															
      arcsec2rad = PI_global / (3600.0D0 * 180.0D0)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! INPUT:
! ----------------------------------------------------------------------
	  
! ----------------------------------------------------------------------
! EOP data
! ----------------------------------------------------------------------
! Select EOP solution
! 1. IERS C04
! 2. IERS RS/PC Daily (finals2000A.daily)
! 3. IGS ultra-rapid ERP + IERS RS/PC Daily (dX,dY)
      EOP_sol = 3
! ----------------------------------------------------------------------
! IERS C04
      if (EOP_sol == 1) then 
      EOP_fname = 'eopc04_08_IAU2000.62-now'
      n_interp = 4   	! number of points for EOP interpolation
	  
! IERS RS/PC Daily
      else if (EOP_sol == 2)  then
      EOP_fname = 'finals2000A.daily'
      n_interp = 4   	! number of points for EOP interpolation
	  
! IGS ultra-rapid ERP
      else if (EOP_sol == 3) then
      ERP_fname = 'igu18543_12.erp'
      EOP_fname = 'finals2000A.daily'
      n_interp = 4   	! number of points for EOP interpolation
      end if
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! IAU Precession-Nutation model:
! 2. IAU2000A:			iau_model = 2000
! 3. IAU2006/2000A:		iau_model = 2006
      iau_model = 2000
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! IERS web calculation service:
! 1. EOC
! 2. RS/PC
      iers_web = 2
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! File names for reading/writing EOM 
      if (iers_web == 1) then
          fname_eom_iers = 'eom_iers1.dat'
      else if (iers_web == 2) then
          fname_eom_iers = 'eom_iers2.dat'
      end if
	  
      fname_eom_iers_ar = 'eom_iers_array.dat'
	  
      if (EOP_sol == 1) then 
          fname_eom_lib = 'eom_eop1.dat'
          fname_dEOM = 'eom1_delta.dat'
      else if (EOP_sol == 2)  then
          fname_eom_lib = 'eom_eop2.dat'
          fname_dEOM = 'eom2_delta.dat'
      else if (EOP_sol == 3) then
          fname_eom_lib = 'eom_eop3.dat'
          fname_dEOM = 'eom3_delta.dat'
      end if 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
      PRINT *,"--------------------------------------------------------"
      PRINT *,"----------------- lib7_trs.f90 -------------------------"
      PRINT *," "
      PRINT *,"--------------------------------------------------------"
      PRINT *,"--------------------- INPUT ----------------------------"
      PRINT *, "EOP data used:", EOP_sol
      PRINT *, "IAU Precession-Nutation model:", iau_model
      PRINT *, "IERS web service center:", iers_web
      PRINT *,"--------------------------------------------------------"
      PRINT *," "
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Numerical Validation test based on IERS web calculation services
! ----------------------------------------------------------------------
! Earth Orientation matrix : ITRF to ICRF
! ----------------------------------------------------------------------
	  
! ----------------------------------------------------------------------
! 1. IERS Earth Orientation Center (EOC) web service
! ----------------------------------------------------------------------
      if (iers_web == 1) then
          time_in = 4 ! TAI 	  
          ! READ fname_eom_iers .dat file
          CALL eom_iers (fname_eom_iers)
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! 2. IERS RS/PC web service
! ----------------------------------------------------------------------
      else if (iers_web == 2) then
          time_in = 3 ! TAI 	  
          ! READ fname_eom_iers .dat file
          CALL eom_iers2 (fname_eom_iers)
! ----------------------------------------------------------------------
      end if

! ----------------------------------------------------------------------
! Write R_iers
      CALL write_array (fname_eom_iers_ar)	  
! Array dimensions
      sz1 = SIZE (wrtArray,DIM=1)
      sz2 = SIZE (wrtArray,DIM=2)
! Allocate R_iers
      ALLOCATE ( R_iers (sz1,sz2) , STAT = AllocateStatus)
      R_iers = wrtArray
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Loop for the time period defined by the IERS web service EOM file
! ----------------------------------------------------------------------
      sz1 = SIZE (R_iers,DIM=1)
      sz2 = SIZE (R_iers,DIM=2)
      Nepochs = sz1 
      ALLOCATE (R_lib(Nepochs,10), STAT = AllocateStatus)  
      ALLOCATE (d_eom(Nepochs,10), STAT = AllocateStatus)  
      
      DO epoch_i = 1 , Nepochs 												
         mjd = R_iers(epoch_i,1)
		 
! ----------------------------------------------------------------------
! "Time Systems" transformation											 
! ----------------------------------------------------------------------
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


! ----------------------------------------------------------------------
! Case 1. EOP by IERS data: EOP reading, interpolation and corrections
! ----------------------------------------------------------------------
! - IERS Earth Orientation Center:			C04 solution
! - IERS Rapid Service/Prediction Center:	finals2000A.daily solution
! ----------------------------------------------------------------------
      IF (EOP_sol == 1 .OR. EOP_sol == 2) THEN  
! EOP corrected
      CALL eop_cor (mjd_TT, EOP_fname, EOP_sol, n_interp, EOP_cr)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Case 2. ERP by IGS ultra-rapid data and dX,dY (IAU 2000A) by IERS RS/PC 
! ----------------------------------------------------------------------
      ELSEIF (EOP_sol == 3)  THEN 
! ERP read and interpolation and dX,dY from finals2000A.daily
      CALL eop_igu (mjd_TT, ERP_fname, EOP_fname, EOP_cr)
! ----------------------------------------------------------------------

      END IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! 1. EOM (1): Initial approach | Earth Orientation Matrix (only)
! ----------------------------------------------------------------------
!      CALL eom (mjd_TT, EOP_fname, n_interp, GCRS2ITRS, ITRS2GCRS )
!      CALL eom2 (mjd_TT, EOP_cr, GCRS2ITRS_2, ITRS2GCRS_2)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 2. ICRF-ITRF transformation matrix (including derivatives)
! ----------------------------------------------------------------------
! CRS-TRS transformation matrix (r,v)
! - CRS2TRS:		GCRS to ITRS transformation matrix (position vector)
! - TRS2CRS:		ITRS to GCRS transformation matrix (position vector)
! - d_CRS2TRS:		Derivative of GCRS to ITRS transformation matrix (velocity vector)
!					v_TRS = TRS2CRS * v_CRS + d_CRS2TRS * r_CRS
! - d_TRS2CRS:		ITRS to GCRS transformation matrix (velocity vector)
!					v_CRS = TRS2CRS * v_TRS + d_TRS2CRS * r_TRS
! ----------------------------------------------------------------------
      CALL crs_trs (mjd_TT, EOP_cr, iau_model, CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      ! R computed by Lib7_trs.f90
      if (iers_web == 1) then
          R_lib(epoch_i,1) = mjd_TAI
      else if (iers_web == 2) then
          R_lib(epoch_i,1) = mjd_UTC
      end if
      !R_lib(epoch_i,1) = mjd
	  
      R_lib(epoch_i,2) = TRS2CRS(1,1)
      R_lib(epoch_i,3) = TRS2CRS(1,2)
      R_lib(epoch_i,4) = TRS2CRS(1,3)
      R_lib(epoch_i,5) = TRS2CRS(2,1)
      R_lib(epoch_i,6) = TRS2CRS(2,2)
      R_lib(epoch_i,7) = TRS2CRS(2,3)
      R_lib(epoch_i,8) = TRS2CRS(3,1)
      R_lib(epoch_i,9) = TRS2CRS(3,2)
      R_lib(epoch_i,10) = TRS2CRS(3,3)
      !PRINT *, R_lib(epoch_i,1)
! ----------------------------------------------------------------------
      END DO 

! ----------------------------------------------------------------------
! WRITE EOM to fname_eom_lib
      ALLOCATE ( wrtArray (Nepochs,10) , STAT = AllocateStatus)
	  wrtArray = R_lib
      CALL write_array (fname_eom_lib)	  
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
	  

! ----------------------------------------------------------------------
! Numerical comparison:  R_lib  Vs.  R_iers
! ----------------------------------------------------------------------
      DO i = 1 , Nepochs
         d_eom(i,1) = R_iers(i,1) - R_lib(i,1)
         d_eom(i,2) = R_iers(i,2) - R_lib(i,2)
         d_eom(i,3) = R_iers(i,3) - R_lib(i,3)
         d_eom(i,4) = R_iers(i,4) - R_lib(i,4)
         d_eom(i,5) = R_iers(i,5) - R_lib(i,5)
         d_eom(i,6) = R_iers(i,6) - R_lib(i,6)
         d_eom(i,7) = R_iers(i,7) - R_lib(i,7)
         d_eom(i,8) = R_iers(i,8) - R_lib(i,8)
         d_eom(i,9) = R_iers(i,9) - R_lib(i,9)
         d_eom(i,10) = R_iers(i,10) - R_lib(i,10)
         !PRINT *, "d_eom", d_eom(i,1:10)		 
      END DO
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Statistics (RMS, min, max) 
! ----------------------------------------------------------------------
! EOM elements individually
      ALLOCATE (array1(Nepochs), STAT = AllocateStatus)
      DO i = 2 , 10
	     array1 = d_eom(:,i)
         CALL stats (stat_q) 
         PRINT *, "d_eom(i)", i	  
         PRINT *, "RMS",  stat_q(1)	  
         PRINT *, "Mean", stat_q(4)	  
         PRINT *, "Max",  stat_q(2)	  
         PRINT *, "Min",  stat_q(3)	  
         PRINT *, " "	  
      END DO
      DEALLOCATE (array1,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
! Overall numerical differences
      ALLOCATE (array1(Nepochs * 9), STAT = AllocateStatus)
      i2 = 0
      DO j = 2 , 10 
	     DO i = 1 , Nepochs
	        i2 = i2 + 1
	        array1(i2) = d_eom(i,j)
	        !array1(i2) = abs (d_eom(i,j))
	        !PRINT *,"array1(i2)",j,i,i2,array1(i2)
         END DO
	  END DO
      CALL stats (stat_q) 
      PRINT *, "Stats(overall)"	  
      PRINT *, "RMS",  stat_q(1)	  
      PRINT *, "Mean", stat_q(4)	  
      PRINT *, "Max",  stat_q(2)	  
      PRINT *, "Min",  stat_q(3)	  
      DEALLOCATE (array1,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! WRITE numerical differences to fname_dEOM
      ALLOCATE (wrtArray(Nepochs,10), STAT = AllocateStatus)
	  wrtArray = d_eom
      CALL write_array (fname_dEOM)	  
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      DEALLOCATE (R_lib,   STAT = DeAllocateStatus)
      DEALLOCATE (R_iers,   STAT = DeAllocateStatus)
      DEALLOCATE (d_eom,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

	  
      end

	  