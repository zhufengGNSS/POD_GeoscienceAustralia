      program lib7_trs


! ----------------------------------------------------------------------
! Program:	lib7_trs.f90
! ----------------------------------------------------------------------
! Purpose:
!  - Time scales transformation
!  - EOP data
!  - Earth Orientation matrix
!  - ICRF-ITRF transformation matrix (direct/inverse & derivatives)
!  - Orbital frame transformation
!  - Keplerian elements
!  - Yaw attitude modelling
! ----------------------------------------------------------------------
! Called subroutines (major):
! - time_TT.f90, time_TAI.f90, time_GPS.f90, time_UTC.f90 : Time scales
! - eop_cor.f90:	IERS EOP data processing and tidal corrections (eop_interp.f90, interp.f)
! - CRS_TRS.f90:	ICRF-ITRF transformation matrix (direct/inverse & derivatives)
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia            December 2015
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Lib8_eop.f90
! ----------------------------------------------------------------------
! Called subroutines (major):
! - eop_cor.f90:	IERS EOP data processing and tidal corrections (eop_interp.f90, interp.f)
! - eop_cor2.f90:		EOP data corrections (ORTHO_EOP.f,..., rg_zont2.f : Zonal Tides)
! - eop_igu.f90 (eop_igu_int.f90) : IGS ultra-rapid ERP data processing
! - eop_rd.f90:			EOP data format reading 
! - eop_interp.f90:		IERS EOP data interpolation and corrections (interp.f: xp,yp,UT1-UTC: ocean tidal and libration effects)
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_planets
      USE mdl_arr
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
      DOUBLE PRECISION arcsec2rad
      INTEGER (KIND = prec_int4) :: time_in, i
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
      !DOUBLE PRECISION GCRS2ITRS(3,3), ITRS2GCRS(3,3)
      !DOUBLE PRECISION GCRS2ITRS_2(3,3), ITRS2GCRS_2(3,3), dTM(3,3)  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: dif_MJD, mjd_iers, TAIsec_iers, TAI_UTC_iers
      !DOUBLE PRECISION dif_MJD, mjd_iers, TAIsec_iers, TAI_UTC_iers
      DOUBLE PRECISION dif_EOM(3,3), R_iers(3,3)
      DOUBLE PRECISION dif_PNM(3,3), PN_iers(3,3)
      DOUBLE PRECISION dif_X, dif_Y, X_iersweb, Y_iersweb
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: x_TRS, y_TRS, z_TRS, r_TRS(3), Vx_TRS, Vy_TRS, Vz_TRS, v_TRS(3)
      REAL (KIND = prec_d) :: r_CRS(3), v_CRS(3), v_CRS_1(3), v_CRS_2(3)
      REAL (KIND = prec_d) :: r_TRS_inv(3), v_TRS_inv(3), v_TRS_1(3), v_TRS_2(3)
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: kepler1(9), kepler2(6), GM
      REAL (KIND = prec_d) :: r_CRS_kepl(3), v_CRS_kepl(3), dr_kepler(3), dv_kepler(3)
      REAL (KIND = prec_d) :: Rrtn(3,3), Rrtn_inv(3,3), r_rtn(3)
! ----------------------------------------------------------------------
!! attitude_yaw.f90
!      INTEGER (KIND = prec_int8) :: PRN, satblk, orbdir
!      REAL (KIND = prec_d) :: beta, Yangle(2), Yaw
!      INTEGER (KIND = prec_int8) :: eclipsf
! yaw_attitude.f90
      REAL (KIND = prec_d) :: beta, Mangle, Yangle(2), Yangle_eBX, Yangle_nom, Mangle_nom
      REAL (KIND = prec_d) :: Yaw
      REAL (KIND = prec_d) :: Mangle_e, Mrate_e, Ynom_e
      REAL (KIND = prec_d), Dimension(3) :: eBX_nom, eBX_ecl, eBX_rtn
      INTEGER (KIND = 4) :: PRN, PRN_sp3, orbdir, satblk, eclipsf
! ----------------------------------------------------------------------
! lib3_planets : Variables declaration
      REAL (KIND = prec_d) :: r_sun_crs(3), r_sun_trs(3)
      DOUBLE PRECISION  JD_TT, R_eph(6)
      INTEGER  NTARG, NCTR
      CHARACTER (LEN=100) :: fname_header, fname_data, fname_out	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: Rcrf_bff(3,3), Rrtn_bff(3,3), R_bf_crf(3,3), e_CRF(3), e_BF(3)
      REAL (KIND = prec_d) :: Rtest_I(3,3)
      INTEGER (KIND = prec_int4) :: AllocateStatus, DeAllocateStatus




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
      time_in = 2
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
! State Vector in ITRF
! ----------------------------------------------------------------------
!*  2015  6  1  0  0  0.00000000
!PG01  17835.809256   4999.483252  19053.904501     -4.877550  7  5  6 122       
!*  2015  6  1  0 15  0.00000000
!PG01  18785.611653   6871.764067  17516.389859     -4.877073  6  5  5 120       

! Position vector
      x_TRS = 17835.809256D0
      y_TRS = 4999.483252D0
      z_TRS = 19053.904501D0
	  r_TRS = (/ x_TRS, y_TRS, z_TRS /)
! Velocity vector
      Vx_TRS = (18785.611653D0 - 17835.809256D0) / (15D0 * 60D0)
      Vy_TRS = (6871.764067D0 - 4999.483252D0) / (15D0 * 60D0)
      Vz_TRS = (17516.389859D0 - 19053.904501D0) / (15D0 * 60D0)
	  v_TRS = (/ Vx_TRS, Vy_TRS, Vz_TRS /)

! Conversion: Km to m, Km/sec to m/s
      r_TRS = r_TRS * 1.0D3
      v_TRS = v_TRS * 1.0D3
      PRINT *, "r_TRS", r_TRS
      PRINT *, "v_TRS", v_TRS	  
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Satellite ID (PRN, Block); Required for the yaw attitude modelling
! ----------------------------------------------------------------------
! SV PRN NUMBER (.LE.32 FOR GPS, .GT.32 FOR GLONASS)
      PRN = 01	
! SV BLOCK:  1=I, 2=II, 3=IIA, IIR=(4, 5) IIF=6
      satblk = 6
! DIRECTION OF PROCESSING (1=FORWARD, -1=BACKWARD)
      orbdir = 1 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! EOP data
! ----------------------------------------------------------------------
! EOP solutions options
! 1. IERS C04
! 2. IERS RS/PC Daily (finals2000A.daily)
! 3. IGS ultra-rapid ERP + IERS RS/PC Daily (dX,dY)
      EOP_sol = 1
! ----------------------------------------------------------------------
! EOP data files name:
!
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
! Gravitational constant
      GM = 0.39860044150D+15
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
      PRINT *,"--------------------------------------------------------"
      PRINT *,"------------------ lib7_trs.f90 ------------------------"
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
! EOP data
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

	  ! Numerical differences: EOP_sol3 - EOP_sol2 
      CALL eop_cor (mjd_TT, EOP_fname, 2, n_interp, EOP_cr2)
      !PRINT *,"delta EOP: IGU-IERS" 
      !PRINT *, EOP_cr - EOP_cr2
! ----------------------------------------------------------------------
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
! Reference Frames transformation
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ITRF to ICRF
! ----------------------------------------------------------------------
! r_CRS = TRS2CRS * r_TRS
      CALL matrix_Rr (TRS2CRS,r_TRS , r_CRS)
	  
! v_CRS = TRS2CRS * v_TRS + d_TRS2CRS * r_TRS
      CALL matrix_Rr (TRS2CRS,v_TRS , v_CRS_1)
      CALL matrix_Rr (d_TRS2CRS,r_TRS , v_CRS_2)
      v_CRS = v_CRS_1 + v_CRS_2
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ICRF to ITRF (Test Inverse transformation) 
! ----------------------------------------------------------------------
! r_TRS_inv
      CALL matrix_Rr (CRS2TRS,r_CRS , r_TRS_inv)
! v_TRS_inv
! v_TRS = TRS2CRS * v_CRS + d_CRS2TRS * r_CRS
      CALL matrix_Rr (CRS2TRS,v_CRS , v_TRS_1)
      CALL matrix_Rr (d_CRS2TRS,r_CRS , v_TRS_2)
      v_TRS_inv = v_TRS_1 + v_TRS_2
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Kepler elements
! ----------------------------------------------------------------------
      CALL kepler_z2k (r_CRS,v_CRS,GM , kepler1)

	  kepler2(1:5) = kepler1(1:5)
	  kepler2(6)   = kepler1(7)
      CALL kepler_k2z (kepler2,GM , r_CRS_kepl,v_CRS_kepl)
      
      ! Inverse test : Numerical differences
      dr_kepler = r_CRS_kepl - r_CRS
      dv_kepler = v_CRS_kepl - v_CRS
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Orbital Frame | RTN: Radial, Along (Tangential), Cross (Normal)
! ----------------------------------------------------------------------
      CALL orb_frame(r_CRS, v_CRS, Rrtn) 		! GCRF to Orbital frame
      CALL matrix_Rr (Rrtn,r_CRS , r_rtn) 		! r

! Orbital frame to GCRF : Inverse matrix
      CALL matrix_inv3 (Rrtn, Rrtn_inv)
      Rtest_I = MATMUL(Rrtn, Rrtn_inv)
	  !PRINT *,"Rrtn_inv",Rrtn_inv
      !PRINT *,"I test: Rrtn",Rtest_I

! Wrong approach via transpose matrix
	  !Rrtn_inv = Transpose(Rrtn)	
      !Rtest_I = MATMUL(Rrtn, Rrtn_inv)
	  !PRINT *,"Rrtn_inv",Rrtn_inv
      !PRINT *,"I test: Rrtn",Rtest_I
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Yaw attitude modelling
! ----------------------------------------------------------------------
      PRINT *,"--------------------------------------------------------"
      PRINT *,"-------------------- Attitude --------------------------"

	  
! ----------------------------------------------------------------------
! Planetary orbit data : lib3_planets.f90
! ----------------------------------------------------------------------
! INPUT (lib3_planets.f90)
! ----------------------------------------------------------------------
      fname_out = 'DE.430' 
      fname_header = 'header.430_229'
      fname_data = 'ascp1950.430'
      CALL CATfile (fname_header,fname_data,fname_out)
      CALL asc2eph (fname_out)
! ----------------------------------------------------------------------
! Center celestial body : Earth
      NCTR = 3 
! Target celestial body : Sun
      NTARG = 11 
! ----------------------------------------------------------------------
! JD at the input epoch (in TT)
      JD_TT = mjd_TT + 2400000.5D0
      PRINT *, "JD_TT:", JD_TT	   
! ----------------------------------------------------------------------
! End of INPUT (lib3_planets.f90)
! ----------------------------------------------------------------------	  

! ----------------------------------------------------------------------
! Sun state vector in ICRF (in KM)
      CALL  PLEPH ( JD_TT, NTARG, NCTR, R_eph )
! KM to M
	  r_sun_crs(1) = R_eph(1) * 1000D0
	  r_sun_crs(2) = R_eph(2) * 1000D0
	  r_sun_crs(3) = R_eph(3) * 1000D0
! ----------------------------------------------------------------------
! Sun state vector in ITRF
      CALL matrix_Rr (CRS2TRS,r_sun_crs , r_sun_trs) 
! ----------------------------------------------------------------------
! End of lib3_planets.f90
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Attitude Yaw angle
! ----------------------------------------------------------------------
      PRINT *,"---------------- yaw_attitude.f90 ----------------------"
	  
! ----------------------------------------------------------------------
! Beta angle
      CALL beta_angle (r_CRS, v_CRS, r_sun_crs, beta)
      PRINT *,"beta (deg)", beta
      PRINT *," " 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Yaw angle: Eclipse model (IECLIPS = 1 or 2) or nominal (IECLIPS = 0) 
! ----------------------------------------------------------------------	  
! ICRF	
! Nominal and Eclipsing Yaw-attitude
      CALL yaw_attitude (mjd_GPS , r_CRS, v_CRS, r_sun_crs, beta, PRN, satblk, orbdir, eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle, Mangle_e, Mrate_e, Ynom_e)

!	  ! Write	  
!      attitude_icrf(epoch_wrt , :) = (/ mjd_GPS, time_write, eclipsf * 1.0D0, beta, Mangle, Yangle(1), Yangle(2) /)
!      eclips_icrf(epoch_wrt , :) = (/ mjd_GPS, time_write, eclipsf * 1.0D0, eBX_ecl(1), eBX_ecl(2), eBX_ecl(3), Mrate_e, Mangle_e, Ynom_e, Yangle(2) /)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Yaw angle (in degrees) based on the eclipse status (nominal or eclipsing)
      If (eclipsf == 0) then
		Yaw = Yangle(1)
      Else If (eclipsf == 1 .OR. eclipsf == 2) then
		Yaw = Yangle(2)
      End If
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Transformation matrix between Inertial frame and Body-Fixed frame
	  CALL crf_bff (r_CRS, v_CRS, Yaw, Rcrf_bff, Rrtn_bff)
! ----------------------------------------------------------------------

! Apply transofrmation from inertial frame to body-fixed frame  	  
      e_CRF = r_CRS * ( 1D0 / sqrt(r_CRS(1)**2 + r_CRS(2)**2 + r_CRS(3)**2) )	
      CALL matrix_Rr (Rcrf_bff, e_CRF, e_BF) 


 
	  
	  
	  
	  

      PRINT *,"--------------------------------------------------------"
      PRINT *,"--------------------------------------------------------"
! ----------------------------------------------------------------------
	  PRINT *, "r_TRS", r_TRS
      PRINT *, "v_TRS", v_TRS

	  PRINT *, "r_CRS", r_CRS
      PRINT *, "v_CRS", v_CRS
! ----------------------------------------------------------------------
      !PRINT *, "CRS2TRS", CRS2TRS
      !PRINT *, "TRS2CRS", TRS2CRS
      !PRINT *, "d_CRS2TRS", d_CRS2TRS
      !PRINT *, "d_TRS2CRS", d_TRS2CRS  
! ----------------------------------------------------------------------  
      PRINT *,"------ Test Inverse transformation (ICRF to ITRF) ------"
      PRINT *, "r_TRS_inv - r_TRS", r_TRS_inv - r_TRS
      PRINT *, "v_TRS_inv - v_TRS", v_TRS_inv - v_TRS
      PRINT *,"--------------------------------------------------------"
! ----------------------------------------------------------------------
	  !PRINT *, "kepler1", kepler1
	  !PRINT *, "kepler2", kepler2
	  !PRINT *, "r_CRS_kepl", r_CRS_kepl
	  !PRINT *, "v_CRS_kepl", v_CRS_kepl
      !PRINT *,"-------- Test Inverse conversion: Kepler to r,v --------"
	  !PRINT *, "dr_kepler", dr_kepler
      !PRINT *, "dv_kepler", dv_kepler
      !PRINT *,"--------------------------------------------------------"
! ----------------------------------------------------------------------
      PRINT *, "Rrtn", Rrtn
      PRINT *, "r_rtn", r_rtn
      PRINT *,"--------------------------------------------------------"
! ----------------------------------------------------------------------
      PRINT *,"Yaw (deg)", Yaw
      PRINT *,"Rcrf_bff", Rcrf_bff 
      PRINT *,"Rrtn_bff", Rrtn_bff 
	  PRINT *, "eBX_nom", eBX_nom
	  PRINT *, "eBX_ecl", eBX_ecl
	  PRINT *, "e_BF", e_BF
      PRINT *,"--------------------------------------------------------"
      PRINT *," "



	  
      end

	  