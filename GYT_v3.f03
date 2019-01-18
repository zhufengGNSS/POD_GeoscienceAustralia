      program GYT


! ----------------------------------------------------------------------
! Program:	GYT.f03
! ----------------------------------------------------------------------
! Purpose:
!  Numerical analysis of the GNSS yaw-attitude modelling and phase wind-up 
!  correction over long-time periods including nominal and eclipse seasons.
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia                July 2016
! ----------------------------------------------------------------------
! Last modified
! - Dr. Thomas Papanikolaou, 21 November 2018:
!	Upgrade from Fortran 90 to Fortran 2003 for calling the modified 
!   subroutines (upgraded to F03) of the POD code package 
! ----------------------------------------------------------------------

      USE mdl_precision
      USE mdl_num
      USE mdl_planets
      !USE mdl_write
      !USE mdl_arr
      !USE mdl_eop
      USE m_eop_data
      USE m_eop_cor
      USE m_eop_igu
      USE m_interporb
      USE m_keplerorb
      USE m_writearray
      USE m_sp3
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int4) :: time_in, i
      INTEGER IY, IM, ID, J_flag
      DOUBLE PRECISION DJM0, sec, FD
      REAL (KIND = prec_d) :: mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC, mjd_UT1
! ----------------------------------------------------------------------
      CHARACTER (LEN=50) :: EOP_fname
      CHARACTER (LEN=50) :: ERP_fname
      INTEGER (KIND = prec_int4) :: EOP_Nint
      INTEGER (KIND = prec_int1) :: EOP_sol
      INTEGER (KIND = prec_int2) :: iau_model
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: EOP_days
      REAL (KIND = prec_d) :: EOP_cr(7)
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: mjd_UTC_day
      REAL (KIND = prec_d) :: mjd_igu
      DOUBLE PRECISION CRS2TRS(3,3), TRS2CRS(3,3), d_CRS2TRS(3,3), d_TRS2CRS(3,3)
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: x_TRS, y_TRS, z_TRS, Vx_TRS, Vy_TRS, Vz_TRS
      REAL (KIND = prec_d) :: r_TRS(3), v_TRS(3)
      REAL (KIND = prec_d) :: r_CRS(3), v_CRS(3), v_CRS_1(3), v_CRS_2(3)
      REAL (KIND = prec_d) :: v_TRS_1(3), v_TRS_2(3)
! ----------------------------------------------------------------------
! GNSS satellite numbers - PRN
      CHARACTER (LEN=1) :: yaw_mod
      CHARACTER (LEN=3) :: PRN
      CHARACTER (LEN=1) :: GNSSid
      INTEGER (KIND = 4) :: PRN_sp3, PRN_eclips
      INTEGER (KIND = 4) :: satbf, satblk
      INTEGER (KIND = 4) :: eclipsf, orbdir
      INTEGER (KIND = prec_int2) :: ios  
      CHARACTER (LEN=100) :: fmt_line
! ----------------------------------------------------------------------
! yaw_attitude.f90
      REAL (KIND = prec_d) :: beta, Mangle, Yangle(2), Yangle_eBX, Yangle_nom, Mangle_nom
      REAL (KIND = prec_d) :: Mangle_e, Mrate_e, Ynom_e
      REAL (KIND = prec_d), Dimension(3) :: eBX_nom, eBX_ecl, eBX_rtn, eBX_nom_ITRF, eBX_ecl_ITRF
! ----------------------------------------------------------------------
! lib3_planets : Variables declaration
      DOUBLE PRECISION  ET, JD_TT, R(6)
      INTEGER  NTARG, NCTR
      CHARACTER (LEN=100) :: fname_header, fname_data, fname_out
      REAL (KIND = prec_d) :: r_sun_crs(3), r_sun_trs(3), v_sun_crs(3)
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: GPS_week, GPSweek_mod1024
      REAL (KIND = prec_d) :: GPS_wsec, GPS_day
      REAL (KIND = prec_d) :: sec_day, time_write
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: data_opt
      INTEGER (KIND = prec_int8) :: interv
      INTEGER (KIND = prec_int8) :: Nepochs, epoch_i, sz1, sz2 
      INTEGER (KIND = prec_int8) :: epoch_wrt 
      CHARACTER (LEN=300) :: fname_orb, fname_write, fname_orbint
      CHARACTER (LEN=300) :: fname_orb_0, fname_orb_1, fname_orb_2
      REAL (KIND = prec_q) :: CPU_t0, CPU_t1
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: attitude_icrf
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: eclips_icrf
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: phasewup_icrf, phasewup_itrf
! ----------------------------------------------------------------------
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: orbint
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: orbkepler
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: orb_Crtn, orb_Trtn
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus, pflag, time_tag
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: Vi(3), r_t1(3), r_t2(3), dt
      REAL (KIND = prec_d) :: Vxyz(3), Vrtn(3), Vxyz_norm, Vrtn_norm, eV_xyz(3), eV_rtn(3)	
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: data_rate, Nepoch_int, Norb_int, NPint, NPint_1, NPint_2, i1, i2, ti, j
	  REAL (KIND = prec_d) :: SecDay_i, MJD_i, MJD_tint, tint, Xint, Xint_1, Yint, Yint_1, Zint, Zint_1, X,Y,Z, Vx,Vy,Vz
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: GM_Earth
      INTEGER (KIND = prec_int4) :: Zo_el
      INTEGER (KIND = prec_int8) :: Ndays
	  REAL (KIND = prec_d) :: Zo(6), Sec0, MJDo, MJDyear
	  DOUBLE PRECISION ERM(3,3)
! ----------------------------------------------------------------------
! yawdyn.f90
      REAL (KIND = prec_d) :: beta0
      INTEGER (KIND = prec_int8) :: dparam
! BDS
      CHARACTER (LEN=5) :: BDSorbtype
      REAL (KIND = prec_d) :: beta_t0 
      INTEGER (KIND = prec_int2) :: BetaP
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Start of INPUT Configuration parameters:
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Yaw attitude model/method
! ----------------------------------------------------------------------	    
! 1. Yaw attitude modelling applied based on the models adopted for each GNSS constellation
!    Set yaw_mod variable to the GNSSid (GNSS constellation id letter)
!    GNSS constellation id letters
!    G: GPS
!    R: GLONASS
!    E: Galileo
!	 C: BDS (BeiDou)
!	 J: QZSS
!	 S: SBAS
! ----------------------------------------------------------------------
! 2. Dynamic Yaw-steering method 
!    Set yaw_mod variable to 'D'
! ----------------------------------------------------------------------
yaw_mod = 'G'
!yaw_mod = 'R'
!yaw_mod = 'E'
!yaw_mod = 'C'
!yaw_mod = 'D'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! GNSS PRN number (Constellation ID + Number)
! ----------------------------------------------------------------------
!PRN = 'G21'  ! IIR
!PRN = 'G30'  ! IIA
PRN = 'G03'  ! IIF
!PRN = 'C13'  ! MEO Dynamic-Yaw
!PRN = 'C06'  ! IGSO 
!PRN = 'C14'   ! MEO
!PRN = 'C12'    
!PRN = 'C09' 
!PRN = 'C11'    
!PRN = 'E11'  ! Galileo 
!PRN = 'E12'  ! Galileo 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------	  
! Satellite Body-fixed frame definition type
! ----------------------------------------------------------------------	  
! satbf = 1 : Body-fixed frame according to the IGS Conventions; Cases: GPS Block II,IIA,IIF, GLONASS, BeiDou  
! satbf = 2 : Body-fixed frame X,Y axes reversal; Cases: Galileo, GPS Block IIR 
satbf = 1
! ----------------------------------------------------------------------	


! ----------------------------------------------------------------------
! GPS case only:
! ----------------------------------------------------------------------
! Satellite Block ID (according to the numbering adopted in eclips.f)
! satblk:	1=I, 2=II, 3=IIA, IIR=(4, 5), IIF=6
satblk = 6
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Beidou case only:
! ----------------------------------------------------------------------
! Orbit type:  MEO, IGSO, GEO
!BDSorbtype = 'IGSO'
BDSorbtype = 'MEO'
! ----------------------------------------------------------------------
! Approach for the beta angle computation at the latest epoch that the orbital angle M is equal to 90 degrees
! BetaP = 1 : Orbit backward prediction based on Keplerian elements  
! BetaP = 2 : Orbit backward computation based on numerical interpolation of the precise orbit data sp3 
BetaP = 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Dynamic yaw-steering method
! ----------------------------------------------------------------------
! yawdyn.f90 subroutine input arguments to be configured
! ----------------------------------------------------------------------	  
! Beta angle condition limit for dyanmic yaw steering method
beta0 = 2.0D0
! Design parameter
dparam = 258
! ----------------------------------------------------------------------	  


! ----------------------------------------------------------------------
! GNSS Orbit data
! ----------------------------------------------------------------------
! 2. Orbit based on Lagrange interpolation of sp3 data files 
data_opt = 2

! ----------------------------------------------------------------------
! Orbit Interpolation parameters
! ----------------------------------------------------------------------
! Orbit interpolation step in sec
interv = 30
! Number of data points used in Lagrange interpolation of GNSS orbit data (sp3)   
NPint = 12
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit sp3 file name
! ----------------------------------------------------------------------
fname_orb = 'igs18885.sp3'	! 18/03/2016  GPS IIF
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
EOP_Nint = 4   	! number of points for EOP interpolation
! ----------------------------------------------------------------------
! EOP data file name:

! IERS C04
if (EOP_sol == 1) then 
	EOP_fname = 'eopc04_14_IAU2000.62-now'   ! 'eopc04_08_IAU2000.62-now'	  
! IERS RS/PC Daily
else if (EOP_sol == 2)  then
	EOP_fname = 'finals2000A.daily'  
! IGS ultra-rapid ERP
else if (EOP_sol == 3) then
	ERP_fname = 'igu18543_12.erp'
	EOP_fname = 'finals2000A.daily'
end if
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! End of Major INPUT Configuration
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Minor (additional) INPUT Configuration
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! IAU Precession-Nutation model:
! 2. IAU2000A:			iau_model = 2000
! 3. IAU2006/2000A:		iau_model = 2006
      iau_model = 2000
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Planetary/Lunar Ephemeris data
! ----------------------------------------------------------------------
      fname_out = 'DE.430' 
      fname_header = 'header.430_229'
      fname_data = 'ascp1950.430'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Time Tag to be written in the out files (2nd collumn)
! 1. Sec of Day
! 2. Hours since 0h of day
! 3. Sec of GPS Week
      time_tag = 1
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
! End of Minor INPUT Configuration
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
      PRINT *,"--------------------------------------------------------"
      PRINT *,"-------------------- GYT v.3 ---------------------------"
      PRINT *," "
      PRINT *,"--------------------- INPUT ----------------------------"
      PRINT *, "Yaw attitude model/method: ", yaw_mod
      PRINT *, "PRN: ", PRN
      PRINT *, "GPS Block", satblk
      PRINT *," "
      PRINT *, "Computations interval", interv
      PRINT *, "GNSS orbit data case:", data_opt
      PRINT *, "GNSS orbit data file: ", TRIM(fname_orb)
      PRINT *," "
      PRINT *, "EOP solution:", EOP_sol
      PRINT *, "EOP data file:", EOP_fname
      PRINT *, "IAU Precession-Nutation model: IAU",iau_model
      PRINT *," "
      PRINT *, "Planetary Ephemeris: ", TRIM (fname_data)
      PRINT *,"--------------------------------------------------------"
      PRINT *," "
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Satellite PRN number 
! ----------------------------------------------------------------------
! Consisted by the GNSS constellation id letter (GNSSid) and the number (PRN_sp3)  
! ----------------------------------------------------------------------
! GNSSid:		GNSS constellation id letter
! 				G: GPS
! 				R: GLONASS
! 				E: Galileo
! 				C: BDS (BeiDou)
! 				J: QZSS
! 				S: SBAS
! ----------------------------------------------------------------------
! PRN_sp3:		PRN numbering as adopted in the sp3 format (Numerical value after the GNSS constellation id letter)
! ----------------------------------------------------------------------
fmt_line = '(A1,I2.2)'
READ (PRN, fmt_line , IOSTAT=ios) GNSSid, PRN_sp3
!print *, "GNSSid ", GNSSid
!print *, "PRN_sp3 ", PRN_sp3


! ----------------------------------------------------------------------
! PRN_eclips:	PRN number as adopted in the eclips.f subroutine numbering
!			 	GPS: 		1-32 
! 				GLONASS:	33-64
!				Galileo:	65-100
!				BeiDou:		101-136
! ----------------------------------------------------------------------
If (GNSSid == 'G') then
    PRN_eclips = PRN_sp3
Else if (GNSSid == 'R') then
    PRN_eclips = PRN_sp3 + 32
Else if (GNSSid == 'E') then
    PRN_eclips = PRN_sp3 + 64
Else if (GNSSid == 'C') then
    PRN_eclips = PRN_sp3 + 100
End If
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! GNSS orbit interpolation 
! ----------------------------------------------------------------------
! Case 2: Orbit based on Lagrange interpolation of sp3 data
!Else if (data_opt == 2) then
Call interp_orb (fname_orb, PRN, interv, NPint, orbint)
! ----------------------------------------------------------------------
!PRINT *,"Orbit: Ok"

	
! ----------------------------------------------------------------------
! Planetary/Lunar Ephemeris (DE): Data processing
! ----------------------------------------------------------------------
      CALL CATfile (fname_header,fname_data,fname_out)
      CALL asc2eph (fname_out)
! Center celestial body : Earth
      NCTR = 3 
! Target celestial body : Sun
      NTARG = 11 
! ----------------------------------------------------------------------
	
	
! ----------------------------------------------------------------------
! EOP data
! ----------------------------------------------------------------------
mjd = orbint(1,1)
Call eop_data (mjd, EOP_fname, EOP_sol, EOP_Nint , EOP_days)
! ----------------------------------------------------------------------
!print *, "eop_data IN ", mjd, EOP_fname, EOP_sol, EOP_Nint 
!print *, "EOP_days ", EOP_days 
	
	
! ----------------------------------------------------------------------
! Loop performing computations per epoch
! ----------------------------------------------------------------------
! CPUT time
      CALL cpu_time (CPU_t0)

     Nepochs = size(orbint, DIM = 1) !- 1								
	 print *,"Nepochs", Nepochs

! Allocate arrays
      ALLOCATE (attitude_icrf( Nepochs, 7), STAT = AllocateStatus)
      ALLOCATE (eclips_icrf(Nepochs, 9), STAT = AllocateStatus)
      ALLOCATE (orb_icrf(Nepochs, 8), STAT = AllocateStatus)
      ALLOCATE (orb_itrf(Nepochs, 8), STAT = AllocateStatus)
      !ALLOCATE (orb_Crtn(Nepochs, 8), STAT = AllocateStatus)
      !ALLOCATE (phasewup_icrf(Nepochs, 5), STAT = AllocateStatus)
      !ALLOCATE (phasewup_itrf(Nepochs, 5), STAT = AllocateStatus)  
	  

! ----------------------------------------------------------------------
! Delete previous records in the output data files before write data 
! ----------------------------------------------------------------------
fname_write = 'attitude.out'
OPEN  (UNIT=7, FILE=fname_write)
CLOSE (UNIT=7, STATUS="DELETE")

fname_write = 'bodyX.out'
OPEN  (UNIT=7, FILE=fname_write)
CLOSE (UNIT=7, STATUS="DELETE")

fname_write = 'orb_icrf.out'
OPEN  (UNIT=7, FILE=fname_write)
CLOSE (UNIT=7, STATUS="DELETE")

fname_write = 'orb_itrf.out'
OPEN  (UNIT=7, FILE=fname_write)
CLOSE (UNIT=7, STATUS="DELETE")
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Epoch-wise computations
! ----------------------------------------------------------------------
epoch_i = 0
DO epoch_i = 1 , Nepochs 
! Sec of day
      sec_day = orbint(epoch_i , 2)  
! MJD (including the fraction of the day)
      mjd = orbint(epoch_i , 1)
  
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
! GPS week and Seconds of GPS week
      CALL time_GPSweek (mjd_GPS , GPS_week, GPS_wsec, GPSweek_mod1024)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Orbit based on precise orbit data (interpolated orbit)
! ----------------------------------------------------------------------
! Satellite state vector (in ITRF)
      r_TRS = orbint(epoch_i , 3:5)
      v_TRS = orbint(epoch_i , 6:8)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! EOP data processing
! ----------------------------------------------------------------------
! Case 1. EOP by IERS data: EOP reading, interpolation and corrections
! ----------------------------------------------------------------------
! - IERS Earth Orientation Center:			C04 solution
! - IERS Rapid Service/Prediction Center:	finals2000A.daily solution
! ----------------------------------------------------------------------
      IF (EOP_sol == 1 .OR. EOP_sol == 2) THEN  
		CALL eop_cor (mjd_TT, EOP_days, EOP_sol, EOP_Nint , EOP_cr)
! ----------------------------------------------------------------------
! Case 2. ERP by IGS ultra-rapid data (ERP reading and interpolation) +
!		  dX,dY w.r.t. IAU2000A by IERS RS/PC (finals2000A.daily) 
! ----------------------------------------------------------------------
      ELSEIF (EOP_sol == 3)  THEN 
		CALL eop_igu (mjd_TT, ERP_fname, EOP_days, EOP_cr)
      END IF
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! Transformation: ITRF to ICRF
! ----------------------------------------------------------------------
! ICRF-ITRF transformation matrix (including derivatives)
      CALL crs_trs (mjd_TT, EOP_cr, iau_model, CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS)
! r_CRS = TRS2CRS * r_TRS
      CALL matrix_Rr (TRS2CRS,r_TRS , r_CRS)
! v_CRS = TRS2CRS * v_TRS + d_TRS2CRS * r_TRS
      CALL matrix_Rr (TRS2CRS,v_TRS , v_CRS_1)
      CALL matrix_Rr (d_TRS2CRS,r_TRS , v_CRS_2)
      v_CRS = v_CRS_1 + v_CRS_2
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Planetary orbit data : lib3_planets.f90
! ----------------------------------------------------------------------
! JD at the input epoch (in TT)
      JD_TT = mjd_TT + 2400000.5D0

! Sun state vector in ICRF (in Km, Km/sec)
      CALL  PLEPH ( JD_TT, NTARG, NCTR, R )

! Km to m
	  r_sun_crs(1) = R(1) * 1000D0
	  r_sun_crs(2) = R(2) * 1000D0
	  r_sun_crs(3) = R(3) * 1000D0

	  v_sun_crs(1) = R(4) * 1000D0
	  v_sun_crs(2) = R(5) * 1000D0
	  v_sun_crs(3) = R(6) * 1000D0
	  
! Sun state vector in ITRF
      CALL matrix_Rr (CRS2TRS,r_sun_crs , r_sun_trs)
! ----------------------------------------------------------------------
! End of lib3_planets.f90
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Beta angle
      CALL beta_angle (r_CRS, v_CRS, r_sun_crs, beta)
! ----------------------------------------------------------------------

! Time tag to be written in the output files (2nd column) | 1st: MJD, 2nd: Time tag (additional)
      If (time_tag == 1) then
		time_write = sec_day
      Else if (time_tag == 2) then
		time_write = sec_day / 3600D0
      Else if (time_tag == 3) then
		time_write = GPS_wsec
      End if

	  
! ----------------------------------------------------------------------
! Yaw angle computation based on the selected attitude model
! ----------------------------------------------------------------------	    
! ICRF	
If (yaw_mod == 'G' .or. yaw_mod == 'R') Then

! GPS and GLONASS yaw-attitude model

! Eclipsing attitude modelling is computed through the modified version of the eclips.f
! orbdir: Direction of processing (1=FORWARD, -1=BACKWARD); eclips.f input argument
orbdir = 1 

! Nominal and Eclipsing Yaw-attitude
CALL yaw_attitude (mjd_GPS , r_CRS, v_CRS, r_sun_crs, beta, PRN_eclips, satblk, orbdir, & 
				   eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle, Mangle_e, Mrate_e, Ynom_e)
!PRINT *,"Yangle", Yangle(1), Yangle(2)
if (eclipsf == 0) then
    Yangle(2) = Yangle(1)
end if

Else if (yaw_mod == 'E') Then

! Galileo attitude law

! Yaw-steering as provided by the European GNSS Service Centre 
CALL yaw_gal (mjd_GPS, r_CRS, v_CRS, r_sun_crs, beta, eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle)
!PRINT *,"Yangle", Yangle(1), Yangle(2)


Else if (yaw_mod == 'C') Then

! BeiDou attitude law
CALL yaw_bds (mjd_GPS, r_CRS, v_CRS, r_sun_crs, BDSorbtype, satbf, BetaP, NPint, & 
			  beta, beta_t0, eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle)

Else if (yaw_mod == 'D') Then

! Dynamic Yaw-steering method
CALL yawdyn (mjd_GPS, r_CRS, v_CRS, r_sun_crs, satbf, beta0, dparam, beta, eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle)

End If
! ----------------------------------------------------------------------	    
!print *,"Eclipse status: ", eclipsf

! ----------------------------------------------------------------------	    
! Write	results per epoch to arrays  
! ----------------------------------------------------------------------	    
attitude_icrf(epoch_i , :) = (/ mjd_GPS, time_write, eclipsf * 1.0D0, beta, Mangle, Yangle(1), Yangle(2) /)

! BeiDou case: Activate only for storing the beta-t0 (m=90) variations | "beta vs beta-t0"
!attitude_icrf(epoch_i , :) = (/ mjd_GPS, beta_t0, eclipsf * 1.0D0, beta, Mangle, Yangle(1), Yangle(2) /)

eclips_icrf(epoch_i , :) = (/ mjd_GPS, time_write, eclipsf * 1.0D0, & 
			eBX_ecl(1), eBX_ecl(2), eBX_ecl(3), eBX_nom(1), eBX_nom(2), eBX_nom(3) /)
! ----------------------------------------------------------------------	  

! End of yaw-attitude computations
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbits to be written in output files
! ----------------------------------------------------------------------
! ICRF 
orb_icrf (epoch_i , :) = (/ mjd_GPS, time_write, r_CRS, v_CRS /)
! ITRF
orb_itrf (epoch_i , :) = (/ mjd_GPS, time_write, r_TRS, v_TRS /)
! ----------------------------------------------------------------------

End Do
PRINT *,"Yaw-attitude computations: Completed"
! End of computations
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write to output files
! ----------------------------------------------------------------------
! Yaw-attitude
fname_write = "attitude.out"
Call writearray (attitude_icrf, fname_write)

! Body-X vector	  
fname_write = "bodyX.out"
Call writearray (eclips_icrf, fname_write)
	  
! Orbits
fname_write = "orb_icrf.out"
Call writearray (orb_icrf, fname_write)
fname_write = "orb_itrf.out"
Call writearray (orb_itrf, fname_write)
! ----------------------------------------------------------------------



CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


End program

