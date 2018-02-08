      program lib8_eop


! ----------------------------------------------------------------------
! Program:	lib8_eop.f90
! ----------------------------------------------------------------------
! Purpose:
!  EOP data preprocessing which includes reading, interpolation and 
!  corrections to Earth rotation due to tidal variations
! ----------------------------------------------------------------------
! Called subroutines (major):
! - eop_cor.f90:	IERS EOP data processing and tidal corrections (eop_interp.f90, interp.f)
! - eop_cor2.f90:		EOP data corrections (ORTHO_EOP.f,..., rg_zont2.f : Zonal Tides)
! - eop_igu.f90 (eop_igu_int.f90) : IGS ultra-rapid ERP data processing
! - eop_rd.f90:			EOP data format reading 
! - eop_interp.f90:		IERS EOP data interpolation and corrections (interp.f: xp,yp,UT1-UTC: ocean tidal and libration effects)
! ----------------------------------------------------------------------
! - time_TT.f90, time_TAI.f90, time_GPS.f90, time_UTC.f90 : Time scales change
! - CRS_TRS.f90:	ICRF-ITRF transformation matrix (direct/inverse & derivatives)
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia            December 2015
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
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
      CHARACTER (LEN=50) :: ERP_fname
      REAL (KIND = prec_d) :: mjd_igu
      DOUBLE PRECISION CRS2TRS(3,3), TRS2CRS(3,3), d_CRS2TRS(3,3), d_TRS2CRS(3,3)
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: dif_MJD, mjd_iers, TAIsec_iers, TAI_UTC_iers
      !DOUBLE PRECISION dif_MJD, mjd_iers, TAIsec_iers, TAI_UTC_iers
      DOUBLE PRECISION dif_EOM(3,3), R_iers(3,3)
      DOUBLE PRECISION dif_PNM(3,3), PN_iers(3,3)
      DOUBLE PRECISION dif_X, dif_Y, X_iersweb, Y_iersweb
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
      arcsec2rad = PI_global / (3600.0D0 * 180.0D0)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! INPUT:
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Time System of input epoch:
! ----------------------------------------------------------------------
! 1. TT
! 2. GPS time
! 3. UTC
! 4. TAI
      time_in = 3
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Date
      IY = 2015
      IM = 7
      ID = 22
      sec = 600.0D0  
! ----------------------------------------------------------------------
! MJD of input date including fraction of the day	  
      CALL iau_CAL2JD ( IY, IM, ID, DJM0, mjd, J_flag )
      mjd = mjd + sec / 86400D0
! ----------------------------------------------------------------------
	  
	  
! ----------------------------------------------------------------------
! EOP data
! ----------------------------------------------------------------------
! EOP solutions options
! 1. IERS C04
! 2. IERS RS/PC Daily (finals2000A.daily)
! 3. IGS ultra-rapid ERP + IERS RS/PC Daily (dX,dY)
      EOP_sol = 3
! ----------------------------------------------------------------------
! EOP data files name:

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
      PRINT *,"--------------------------------------------------------"
      PRINT *,"----------------- lib8_trs.f90 -------------------------"
      PRINT *," "
      PRINT *,"--------------------------------------------------------"
      PRINT *,"--------------------- INPUT ----------------------------"
      PRINT *, "mjd", mjd 
      PRINT *, "EOP data:", EOP_sol
      PRINT *, "IAU Precession-Nutation model: IAU",iau_model
      PRINT *,"--------------------------------------------------------"
      PRINT *," "
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------



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
! EOP subroutines
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
      CALL eop_cor2 (mjd_TT , delta_EOP )
      PRINT *, "EOP corrections: delta_EOP based on eop_cor2.f90"
      PRINT *, delta_EOP
      PRINT *," "
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Case 1. EOP by IERS data: EOP reading, interpolation and corrections
! ----------------------------------------------------------------------
! - IERS Earth Orientation Center:			C04 solution
! - IERS Rapid Service/Prediction Center:	finals2000A.daily solution
! ----------------------------------------------------------------------
      IF (EOP_sol == 1 .OR. EOP_sol == 2) THEN  
      CALL eop_cor (mjd_TT, EOP_fname, EOP_sol, n_interp, EOP_cr)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Case 2. ERP by IGS ultra-rapid data (ERP reading and interpolation) +
!		  dX,dY w.r.t. IAU2000A by IERS RS/PC (finals2000A.daily) 
! ----------------------------------------------------------------------
      ELSEIF (EOP_sol == 3)  THEN 
      CALL eop_igu (mjd_TT, ERP_fname, EOP_fname, EOP_cr)
! ----------------------------------------------------------------------
	  ! Numerical differences: EOP_sol3 - EOP_sol2 
      CALL eop_cor (mjd_TT, EOP_fname, 2, n_interp, EOP_cr2)
      !PRINT *,"delta EOP: IGU-IERS" 
      !PRINT *, EOP_cr - EOP_cr2
	  
      END IF
! ----------------------------------------------------------------------
 


! ----------------------------------------------------------------------
! 1. Earth Orientation Matrix (only)
! ----------------------------------------------------------------------
!      CALL eom (mjd_TT, EOP_fname, n_interp, GCRS2ITRS, ITRS2GCRS )
!      CALL eom2 (mjd_TT, EOP_cr, GCRS2ITRS_2, ITRS2GCRS_2)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 2. ICRF-ITRF transformation matrix (including derivatives)
! ----------------------------------------------------------------------
! CRS-TRS transformation matrix for position & velocity vectors (r,v)
! - CRS2TRS:		GCRS to ITRS transformation matrix (position vector)
! - TRS2CRS:		ITRS to GCRS transformation matrix (position vector)
! - d_CRS2TRS:		Derivative of GCRS to ITRS transformation matrix (velocity vector)
!					v_TRS = CRS2TRS * v_CRS + d_CRS2TRS * r_CRS
! - d_TRS2CRS:		ITRS to GCRS transformation matrix (velocity vector)
!					v_CRS = TRS2CRS * v_TRS + d_TRS2CRS * r_TRS
! ----------------------------------------------------------------------
      CALL crs_trs (mjd_TT, EOP_cr, iau_model, CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Numerical test
! ----------------------------------------------------------------------
! Earth Orientation matrix : ITRF to ICRF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1. IERS Earth Orientation Center (EOC) web service
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 2015 07 22  0h 10min 0sec
! ----------------------------------------------------------------------
!        TAI          //          Civil date UTC               // Matrix 
!Modified Julian Date //  an  / mois / jour / hour / min / sec / R11 / R12 / R13 / R21 / R22 / R23 / R31 / R32 / R33 
!57225.007372685184 2015  7 22  0   10 -0.0000001  0.526055144627  0.850449107614  0.001516631053 -0.850450047347  0.526055809555 -0.000046904387 -0.000837722370 -0.001265144657  0.999998848814
!57225.049039351848 2015  7 22  1   10 -0.0000003  0.728733695213  0.684795519017  0.001516111928 -0.684796253795  0.728734580362 -0.000046624912 -0.001136771720 -0.001004250624  0.999998849615
! ----------------------------------------------------------------------
      mjd_iers = 57225.007372685184D0
      R_iers(1,1) =  0.526055144627D0
      R_iers(1,2) =  0.850449107614D0
      R_iers(1,3) =  0.001516631053D0
      R_iers(2,1) = -0.850450047347D0
      R_iers(2,2) =  0.526055809555D0
      R_iers(2,3) = -0.000046904387D0
      R_iers(3,1) = -0.000837722370D0
      R_iers(3,2) = -0.001265144657D0
      R_iers(3,3) =  0.999998848814D0
! ----------------------------------------------------------------------
 

! ----------------------------------------------------------------------
! 2. IERS RS/PC web service
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Numerical comparison	  
! ----------------------------------------------------------------------
      PRINT *,"--- Comparison: Lib7_trs.f90 vs IERS web service ---"
! ----------------------------------------------------------------------
! MJD TAI differences: 1 sec (???)
      dif_MJD = mjd_iers - mjd_TAI !- 1.0D0 / 86400D0
      PRINT *,"mjd_iers - mjd_TAI", dif_MJD
! ----------------------------------------------------------------------
! Earth Orientation matrix: ITRF to ICRF
      PRINT *,"--------------------------------------------------------"
      dif_EOM = R_iers - TRS2CRS
      PRINT *,"EOM differences" 
      PRINT *, dif_EOM 
      PRINT *,"--------------------------------------------------------"
! ----------------------------------------------------------------------
      PRINT *," "


	  
      end

	  