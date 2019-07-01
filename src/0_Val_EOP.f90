      program Val_EOP


! ----------------------------------------------------------------------
! Program:	Val_EOP.f90
! ----------------------------------------------------------------------
! Purpose:
!  Library tasks:
!  - EOP data processing
!  - Earth Orientation matrix computing
!  - ICRF-ITRF transformation matrix (direct/inverse & derivatives)
! ----------------------------------------------------------------------
! Called subroutines (basic):
! - time_TT.f90, time_TAI.f90, time_GPS.f90, time_UTC.f90 : Time scales change
! - eop_cor.f90:	IERS EOP data processing and tidal corrections (eop_interp.f90, interp.f)
! - eop_urapid.f90 (eop_igu_int.f90) : IGS ultra-rapid ERP data processing
! - CRS_TRS.f90:	ICRF-ITRF transformation matrix (direct/inverse & derivatives)
!
! Subroutines (intermediate):
! - eop_rd.f90:			EOP data format reading 
! - eop_interp.f90:		IERS EOP data interpolation and corrections (interp.f: xp,yp,UT1-UTC: ocean tidal and libration effects)
! - eop_cor2.f90:		EOP data corrections (ORTHO_EOP.f, ..., rg_zont2.f : Zonal Tides)
! - eom.f90, eom2.f90:	Earth Orientation Matrix (initial subroutines : Deactivated)
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia               April 2016
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_write
      IMPLICIT NONE

! ----------------------------------------------------------------------
      DOUBLE PRECISION arcsec2rad
      INTEGER (KIND = prec_int4) :: time_in
      INTEGER IY, IM, ID, J_flag
      DOUBLE PRECISION DJM0, sec, FD
      REAL (KIND = prec_d) :: mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC, mjd_UT1, mjd_0
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
      INTEGER (KIND = prec_int4) :: epoch_i, t_step, t_period, Nepochs, sz1, sz2, i, iers_web
	  CHARACTER (LEN=150) :: fname_eom_iers, fname_eom_iers_ar, fname_eom_lib, fname_dEOM 
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus




! ----------------------------------------------------------------------
! PI precision quad															
      arcsec2rad = PI_global / (3600.0D0 * 180.0D0)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! INPUT:
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Time System of input epoch:
! 1. TT
! 2. GPS time
! 3. UTC
! 4. TAI
      time_in = 3
! ----------------------------------------------------------------------
! Date
      IY = 2015
      IM = 7
      ID = 22
      sec = 0.0D0  
! ----------------------------------------------------------------------
! Step interval (in sec)
      t_step = 300
! Time period in sec (e.g. 6 hours)
      t_period = 1 * 6 * 60 * 60
! Computation epochs
      Nepochs = t_period / t_step 
! ----------------------------------------------------------------------
      ALLOCATE (R_lib(Nepochs,10), STAT = AllocateStatus)  
      ALLOCATE (d_eom(Nepochs,10), STAT = AllocateStatus)  
 
 
	  
	  
! ----------------------------------------------------------------------
! EOP data
! ----------------------------------------------------------------------
! Select EOP solution
! 1. IERS C04
! 2. IERS RS/PC Daily (finals2000A.daily)
! 3. IGS ultra-rapid ERP + IERS RS/PC Daily (dX,dY)
      EOP_sol = 1
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
      iers_web = 1
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
      CALL iau_CAL2JD ( IY, IM, ID, DJM0, mjd_0, J_flag )
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Loop for the time period
! ----------------------------------------------------------------------
      DO epoch_i = 1 , Nepochs 												! 1 : t_step : t_period
	     sec = sec + t_step * (epoch_i - 1)
		 
         mjd = mjd_0 + (t_step * epoch_i) / 86400D0 
         !mjd = mjd_0 + (t_step * epoch_i) / 86400D0 + 1D0 / 86400D0 ! Exclude leap second
      !mjd = mjd_0 + sec / 86400D0
      !mjd = mjd + sec / 86400D0 + 1D0 / 86400D0 ! Exclude leap second
! ----------------------------------------------------------------------
       PRINT *, "mjd", mjd

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
! Individual calls to subroutines
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! EOP reading
! ----------------------------------------------------------------------
!      mjd_UTC_day = INT (mjd_UTC)
!      CALL eop_rd (EOP_fname, EOP_sol, mjd_UTC_day , EOP_data)
      !CALL eop_finals2000A (EOP_fname,mjd_UTC_day , EOP_data)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Tidal variations in Earth Rotation due to ocean tides, libration effects, zonal tides
! ----------------------------------------------------------------------
!      CALL eop_cor2 (mjd_TT , delta_EOP )
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Case 1. EOP by IERS data: EOP reading, interpolation and corrections
! ----------------------------------------------------------------------
! - IERS Earth Orientation Center:			C04 solution
! - IERS Rapid Service/Prediction Center:	finals2000A.daily solution
! ----------------------------------------------------------------------
      IF (EOP_sol == 1 .OR. EOP_sol == 2) THEN  
! ----------------------------------------------------------------------
! EOP corrected
      CALL eop_cor (mjd_TT, EOP_fname, EOP_sol, n_interp, EOP_cr)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Case 2. ERP by IGS ultra-rapid data: EOP reading and interpolation
! ----------------------------------------------------------------------
      ELSEIF (EOP_sol == 3)  THEN 
! ----------------------------------------------------------------------
! EOP ultra-rapid: ERP by IGS and dX,dY (IAU 2000A) by IERS RS/PC (finals2000A.daily)
      CALL eop_igu (mjd_TT, ERP_fname, EOP_fname, EOP_cr)
! ----------------------------------------------------------------------
! Numerical differences
      CALL eop_cor (mjd_TT, EOP_fname, 2, n_interp, EOP_cr2)
      PRINT *,"delta EOP: IGU-IERS" 
      PRINT *, EOP_cr - EOP_cr2
! ----------------------------------------------------------------------
      END IF
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! 1. EOM (1): Initial approach | Earth Orientation Matrix (only)
! ----------------------------------------------------------------------
! Earth Orientation Matrix
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
      PRINT *, fname_eom_lib

! ----------------------------------------------------------------------
      ALLOCATE ( wrtArray (Nepochs,10) , STAT = AllocateStatus)
!	  wrtArray = R_lib
!      CALL write_array (fname_eom_lib)	  
!      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
! Array dimensions
!      sz1 = SIZE (R_lib,DIM=1)
!      sz2 = SIZE (R_lib,DIM=2)
!      PRINT *, "sz1,sz2,Nepochs", sz1, sz2, Nepochs
! ----------------------------------------------------------------------
!      ALLOCATE ( wrtArray (sz1,sz2) , STAT = AllocateStatus)
!      ALLOCATE ( wrtArray (Nepochs,10) , STAT = AllocateStatus)
! ----------------------------------------------------------------------
	  wrtArray = R_lib
      CALL write_array (fname_eom_lib)	  
      !CALL write_array (TRIM (fname_eom_lib) )	  
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
	  
	  

! ----------------------------------------------------------------------
! Numerical Validation test
! ----------------------------------------------------------------------
! Earth Orientation matrix : ITRF to ICRF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1. IERS Earth Orientation Center (EOC) web service
! ----------------------------------------------------------------------
      if (iers_web == 1) then
          ! READ fname_eom_iers .dat file
          CALL eom_iers (fname_eom_iers)
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! 2. IERS RS/PC web service
! ----------------------------------------------------------------------
      else if (iers_web == 2) then
          ! READ fname_eom_iers .dat file
          CALL eom_iers2 (fname_eom_iers)
! ----------------------------------------------------------------------
      end if

! ----------------------------------------------------------------------
! Write R_iers
      CALL write_array (fname_eom_iers_ar)	 
      !CALL write_array (TRIM(fname_eom_iers_ar))	 
! Array dimensions
      sz1 = SIZE (wrtArray,DIM=1)
      sz2 = SIZE (wrtArray,DIM=2)
! Allocate R_iers
      ALLOCATE ( R_iers (sz1,sz2) , STAT = AllocateStatus)
      R_iers = wrtArray
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Numerical comparison:  R_lib  Vs.  R_iers
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

         PRINT *, "d_eom", d_eom(i,1:10)		 
		 
      END DO
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Statistics (RMS, min, max) : Overall and each of the EOM elements individually 

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! WRITE numerical differences to fname_dEOM
      ALLOCATE (wrtArray(Nepochs,10), STAT = AllocateStatus)
	  wrtArray = d_eom
      CALL write_array (fname_dEOM)	  
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      !DEALLOCATE (R_lib,   STAT = DeAllocateStatus)
      !DEALLOCATE (R_iers,   STAT = DeAllocateStatus)
      !DEALLOCATE (d_eom,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

	  
      end

	  