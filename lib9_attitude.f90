      program lib9_attitude


! ----------------------------------------------------------------------
! Program:	lib9_attitude.f90
! ----------------------------------------------------------------------
! Purpose:
!  The yaw-attitude modelling is based on models introduced by Kouba (2009),
!  Dilssner (2010) and Dilssner et al. (2011).
!  The subroutine eclips.f, written by Kouba, implements these models and 
!  is used here.
! ----------------------------------------------------------------------
! Called subroutines (major):
! - eclips.f
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia                 May 2016
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_planets
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
      REAL (KIND = prec_q) :: kepler(6), GM
      REAL (KIND = prec_d) :: r_CRS_kepl(3), v_CRS_kepl(3), dr_kepler(3), dv_kepler(3)
      REAL (KIND = prec_d) :: Rrtn(3,3), r_rtn(3), dr_rtn(3), e_r_rtn(3)
      REAL (KIND = prec_d) :: r_sun_crs(3), r_sun_trs(3)
! ----------------------------------------------------------------------
! attitude_yaw.f90
      REAL (KIND = prec_d) :: beta, Yangle(2)
      REAL (KIND = prec_d) :: e_BX(3)
      !INTEGER (KIND = prec_int8) :: PRN, orbdir
      !INTEGER (KIND = prec_int4) :: satblk
      !NTEGER (KIND = prec_int8) :: eclipsf
      INTEGER (KIND = 4) :: PRN, orbdir, satblk, eclipsf
! ----------------------------------------------------------------------
! lib3_planets : Variables declaration
      DOUBLE PRECISION  ET, JD_TT, R(6)
      DOUBLE PRECISION  R_Horizon(6), dR(6)
      INTEGER  NTARG, NCTR
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: GPS_week, GPSweek_mod1024
      REAL (KIND = prec_d) :: GPS_wsec
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: V_dif_test(3), v2_TRS_test(3), v_test(3)

	  
	  

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
      IY = 2010
      IM = 12
      ID = 25
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! - r_TRS:	Satellite position vector (in ITRF)
! - v_TRS:	Satellite velocity vector (in ITRF)
! - sec:	Seconds of Day
! - PRN : 	SV PRN Number (.LE.32 FOR GPS, .GT.32 FOR GLONASS)
! - satblk:	SV BLOCK:  1=I, 2=II, 3=IIA, IIR=(4, 5), IIF=6
! - orbdir:	Direction of processing (1=FORWARD, -1=BACKWARD)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Case 1.
! ----------------------------------------------------------------------
!*  2010 12 25  4 15  0.00000000                                                 
!PG30   7364.279409 -24470.133488   6859.353833    307.722323                    
!*  2010 12 25  4 30  0.00000000                                                 
!PG30   7543.289649 -23467.631104   9499.003931    307.726893                    
! ----------------------------------------------------------------------
      !r_TRS = (/ 7543.289649D0, -23467.631104D0,   9499.003931D0  /)
      !v_TRS = ( (/ 7543.289649D0, -23467.631104D0,   9499.003931D0 /) - (/ 7364.279409D0, -24470.133488D0,   6859.353833D0 /) ) / (15D0 * 60D0)
      !sec = 4D0 * 3600D0 + 30D0 * 60D0 

	  !v2_TRS_test = -1.0D0 * ( (/ 7543.289649D0, -23467.631104D0,   9499.003931D0 /) - (/ 7773.013007D0, -22181.989098D0,  11972.470848D0  /) ) / (15D0 * 60D0) 
	  !V_dif_test = v2_TRS_test - v_TRS
      !PRINT *, "v_TRS", v_TRS
      !PRINT *, "v2_TRS_test", v2_TRS_test
      !PRINT *, "V_dif_test", V_dif_test
	  
! ----------------------------------------------------------------------
!*  2010 12 25  4 45  0.00000000                                                 
!PG30   7773.013007 -22181.989098  11972.470848    307.729435  
!*  2010 12 25  5  0  0.00000000                                                 
!PG30   8085.454927 -20637.254174  14235.574040    307.733014                    
      !r_TRS = (/ 8085.454927D0, -20637.254174D0,  14235.574040D0  /)
      !v_TRS = ( (/ 8085.454927D0, -20637.254174D0,  14235.574040D0  /) - (/ 7773.013007D0, -22181.989098D0,  11972.470848D0  /) ) / (15D0 * 60D0)
      !sec = 5D0 * 3600D0
	  
!*  2010 12 25  5  0  0.0000
!P 30   8085.454948 -20637.254167  14235.574045    307.738186
!V 30   4037.221569  18492.158963  23822.311443 999999.999999
      r_TRS = (/ 8085.454948D0, -20637.254167D0, 14235.574045D0 /)
      v_TRS = (/ 4037.221569D0, 18492.158963D0, 23822.311443D0 /) / 1.0D4
      sec = 5D0 * 3600D0
! ----------------------------------------------------------------------
      PRN = 30	
      satblk = 3
      orbdir = 1 
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Case 2.
! ----------------------------------------------------------------------
!*  2010 12 25  5 30  0.00000000                                                 
!PG16   2131.569026 -24889.697131   8622.500556   -137.196049 
!*  2010 12 25  5 45  0.00000000                                                 
!PG16   2546.871941 -23810.290114  11208.187125   -137.199089                    
      !r_TRS = (/ 2546.871941D0, -23810.290114D0,  11208.187125D0 /)
      !v_TRS = ( (/ 2546.871941D0, -23810.290114D0,  11208.187125D0 /) - (/ 2131.569026D0, -24889.697131D0,   8622.500556D0 /) ) / (15D0 * 60D0)
      !sec = 5D0 * 3600D0 + 45D0 * 60D0
! ----------------------------------------------------------------------
!*  2010 12 25  5 45  0.0000
!P 16   2546.871955 -23810.290200  11208.187136   -137.194149
!V 16   5258.488671  13447.239459  27723.642297 999999.999999
      !r_TRS = (/ 2546.871955D0, -23810.290200D0,  11208.187136D0 /)
      !v_TRS = (/ 5258.488671D0,  13447.239459D0,  27723.642297D0 /) / 1.0D4
      !sec = 5D0 * 3600D0 + 45D0 * 60D0
	  
      !PRN = 16	
      !satblk = 4
      !orbdir = 1 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Case 3.
! ----------------------------------------------------------------------
!*  2010 12 25  2 30  0.00000000                                                 
!PG25  18838.196203 -16068.018598   9602.748162    -59.793733  
!*  2010 12 25  2 45  0.00000000                                                 
!PG25  18364.799222 -14898.816479  12080.717249    -59.795769   
      !r_TRS =   (/ 18364.799222D0, -14898.816479D0,  12080.717249D0 /)
      !v_TRS = ( (/ 18364.799222D0, -14898.816479D0,  12080.717249D0 /) - (/ 18838.196203D0, -16068.018598D0, 9602.748162D0 /) ) / (15D0 * 60D0)
      !sec = 2D0 * 3600D0 + 45D0 * 60D0
! ----------------------------------------------------------------------
!*  2010 12 25  2 45  0.0000
!P 25  18364.799320 -14898.816560  12080.717339    -59.793556
!V 25  -5706.788504  14447.018572  26453.305716 999999.999999
      !r_TRS = (/ 18364.799320D0, -14898.816560D0, 12080.717339D0 /)
      !v_TRS = (/ -5706.788504D0,  14447.018572D0, 26453.305716D0 /) / 1.0D4
      !sec = 2D0 * 3600D0 + 45D0 * 60D0
	  
      !PRN = 25	
      !satblk = 6
      !orbdir = 1 
! ----------------------------------------------------------------------

      !v_test = v_TRS
      !PRINT *,"v_test",v_test
      !PRINT *,"V norm",sqrt( v_test(1)**2 + v_test(2)**2 + v_test(3)**2 )


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
      PRINT *,"---------------- lib9_attitude.f90 ---------------------"
      PRINT *," "
      PRINT *,"--------------------- INPUT ----------------------------"
      PRINT *, "mjd", mjd 
      PRINT *, "IPRN", PRN
      PRINT *, "IBLK", satblk
      PRINT *, "IDIR", orbdir
      PRINT *, "r_TRS", r_TRS
      PRINT *, "v_TRS", v_TRS
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

      CALL time_GPSweek (mjd_GPS , GPS_week, GPS_wsec, GPSweek_mod1024)
      PRINT *, "GPS_week", GPS_week
      PRINT *, "GPS_wsec", GPS_wsec
      PRINT *, " "



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
      END IF
! ----------------------------------------------------------------------
 


! ----------------------------------------------------------------------
! ICRF-ITRF transformation matrix (including derivatives)
! ----------------------------------------------------------------------
      CALL crs_trs (mjd_TT, EOP_cr, iau_model, CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS)
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
      PRINT *,"r_CRS",r_CRS
      PRINT *,"v_CRS",v_CRS


! ----------------------------------------------------------------------
! Orbital Frame | RTN: Radial, Along (Tangential), Cross (Normal)
! ----------------------------------------------------------------------
      CALL orb_frame(r_CRS, v_CRS, Rrtn) 		! in GCRF
      CALL matrix_Rr (Rrtn,r_CRS , r_rtn) 		! r
	  ! Unit vector
      e_r_rtn = ( 1.0D0 / sqrt(r_rtn(1)**2 + r_rtn(2)**2 + r_rtn(3)**2) ) * r_rtn
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Attitude: Yaw attitude of GPS
! ----------------------------------------------------------------------
      PRINT *,"-------------------- Attitude --------------------------"
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Planetary orbit data : lib3_planets.f90
! ----------------------------------------------------------------------

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
! ----------------------------------------------------------------------
! End of INPUT (lib3_planets.f90)
! ----------------------------------------------------------------------	  

! ----------------------------------------------------------------------
! Sun state vector in ICRF (in KM)
      CALL  PLEPH ( JD_TT, NTARG, NCTR, R )
! KM to M
	  r_sun_crs(1) = R(1) * 1000D0
	  r_sun_crs(2) = R(2) * 1000D0
	  r_sun_crs(3) = R(3) * 1000D0
! ----------------------------------------------------------------------
! Sun state vector in ITRF
      CALL matrix_Rr (CRS2TRS,r_sun_crs , r_sun_trs) 
! ----------------------------------------------------------------------
! End of lib3_planets.f90
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! attitude_yaw0.f90
! ----------------------------------------------------------------------
      !PRINT *,"----------------- attitude_yaw.f90 ---------------------"
      !CALL attitude_yaw (mjd_GPS,r_TRS,v_TRS,r_sun_trs, v_CRS)			! in ITRF
      !PRINT *," "
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! attitude_yaw2.f90
! ----------------------------------------------------------------------
      PRINT *,"---------------- attitude_yaw.f90 ----------------------"
	  
! ----------------------------------------------------------------------
! Beta angle
      CALL beta_angle (r_CRS, v_CRS, r_sun_crs, beta)
      PRINT *,"beta (deg)", beta
      PRINT *," " 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Attitude Yaw angle
! ----------------------------------------------------------------------
! ICRF																	| Agreement in signs and close approximation
      PRINT *,"--- ICRF ---"
      CALL attitude_yaw (mjd_GPS , r_CRS, v_CRS, r_sun_crs, beta, PRN, satblk, orbdir, eclipsf, Yangle)
      PRINT *,"Yangle", Yangle
      PRINT *," "

! ITRF																	| NOT (Lead to wrong values of MURATE)
!      PRINT *,"--- ITRF ---"
!      CALL attitude_yaw (mjd_GPS , r_TRS, v_TRS, r_sun_trs, beta, PRN, satblk, orbdir, eclipsf, Yangle)
!      PRINT *,"Yangle", Yangle
!      PRINT *," "

! (Xsat,Xsun) in ITRF and (eVsat) in ICRF according to Kouba (2009)
      PRINT *,"--- r: ITRF | v: ICRF ---"
      CALL attitude_yaw (mjd_GPS , r_TRS, v_CRS, r_sun_trs, beta, PRN, satblk, orbdir, eclipsf, Yangle)
      PRINT *,"Yangle", Yangle
! ----------------------------------------------------------------------
      PRINT *, " "
      PRINT *, "eclipsf", eclipsf
      PRINT *,"--------------------------------------------------------"


! ----------------------------------------------------------------------
!	  sz1 = SIZE (GNSS_id,DIM=1)
!      sz2 = SIZE (GNSS_id,DIM=2)
!	  DO do_i = 1 , sz1
!	     if (GNSS_id(do_i,1) == PRN) then
!		     IBLK(PRN) = GNSS_id(do_i,2)
!         end if
!      End do
! ----------------------------------------------------------------------



	  
      end

	  