MODULE m_pod_gnss


! ----------------------------------------------------------------------
! MODULE: m_pod_gnss.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for Precise Orbit Determination
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	20 May 2019
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE pod_gnss (EQMfname, VEQfname, PRNmatrix, orbits_partials_icrf, orbits_partials_itrf, &
                     orbits_ics_icrf, orbit_resR, orbit_resT, orbit_resN, orbdiff2)

! ----------------------------------------------------------------------
! SUBROUTINE:	pod_gnss.f03
! ----------------------------------------------------------------------
! Purpose:
!  Precise Orbit Determination (POD) of GNSS constellations 
! ----------------------------------------------------------------------
! Input arguments:
! - EQMfname: 	Input configuration file name for the orbit integration of the Equation of Motion  
! - VEQfname: 	Input configuration file name for the orbit integration of the Variational Equations
!
! Output arguments:
! - PRNmatrix:				PRN numbers array e.g. G01, .., G32, E01, .., E30
! - orbits_partials_icrf: 	Satellite Orbits and Partial derivatives of the estimated parameters in inertial frame (ICRF) per satellite per epoch:
!   Format:
! 				Row 1   :: Satellite 1, Epoch 1 :: Format:
!               - Modified Julian Day number (including the fraction of the day) 
!				- Seconds since 00h 
!				- Position vector (m)
!				- Velocity vector (m/sec)
! 				- Partial Derivatives
! 				...
! 				Row N   :: Satellite N, Epoch 1 :: Format as per above
!				Row N+1 :: Satellite 1, Epoch 2 :: Format as per above
!				...
!				Row Nsat*Nepochs :: Satellite N, Epoch Final :: Format as per above 
! - orbits_partials_itrf:   
! - orbits_residuals: 	
!
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	20 May 2019
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_config
      USE m_orbitmain
      USE m_writeorbit_multi
      USE m_orbdet
      USE m_orbext
      USE m_writearray
      USE m_writeorbit
	  USE mdl_planets
	  USE mdl_tides
	  USE mdl_eop
	  USE m_sp3_PRN
	  USE m_write_orb2sp3
      USE m_orbitIC
      USE m_read_satsnx 
	  
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN)  :: EQMfname, VEQfname				
! ----------------------------------------------------------------------
! OUT
	  CHARACTER (LEN=3), ALLOCATABLE, INTENT(OUT) :: PRNmatrix(:)
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: orbits_partials_icrf  
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: orbits_partials_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orbits_ics_icrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orbit_resR  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orbit_resT  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orbit_resN  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: PODfname, ORBMODfname				
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Vres  
      REAL (KIND = prec_d), DIMENSION(3) :: Vrms 	    
	  !REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT
! ----------------------------------------------------------------------
      CHARACTER (LEN=2) :: GNSS_id
	  INTEGER (KIND = prec_int2) :: ORB_mode
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Nsat, isat
	  INTEGER (KIND = prec_int8) :: iepoch, iparam
	  INTEGER (KIND = prec_int8) :: i
	  INTEGER (KIND = prec_int8) :: sz1, sz2, Nepochs, N2_orb, N2_veqSmatrix, N2_veqPmatrix, N2sum , N2ics
      !REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_icrf  
      !REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_itrf  
	  !CHARACTER (LEN=3), ALLOCATABLE :: PRNmatrix(:)
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus  
	  CHARACTER (LEN=3) :: PRN_isat
	  INTEGER :: ios,ios_key
      CHARACTER (LEN=100) :: orbits_fname				
      CHARACTER (LEN=100) :: fname_write				
      CHARACTER (LEN=100) :: filename				
      CHARACTER (LEN=300) :: fname_sp3, ORBpseudobs_fname, ORBEXT_fname				
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: mjd, mjd0, jd0
      REAL (KIND = prec_d) :: Sec_00 	    
	  INTEGER :: year, month, day
	  INTEGER :: Iyear, Imonth, Iday
      INTEGER J_flag
      DOUBLE PRECISION FD
! ----------------------------------------------------------------------
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
      REAL (KIND = prec_d) :: Zo(6) 
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: ORB2sp3_fname				
      INTEGER (KIND = prec_int2) :: sat_vel	  	  
! ----------------------------------------------------------------------
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resR  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resT  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resN  
! ----------------------------------------------------------------------
      !INTEGER (KIND = prec_int8) :: GPS_week, GPSweek_mod1024
      !REAL (KIND = prec_d) :: GPS_wsec, GPS_day
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Ncommon  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_RTN, dorb_Kepler
! ----------------------------------------------------------------------
      !INTEGER (KIND = prec_int2) :: POD_MODE_glb	  	  
      !REAL (KIND = prec_d) :: ORBPRED_ARC_glb
      !REAL (KIND = prec_d) :: orbarc_sum
      !INTEGER (KIND = prec_int2) :: IC_MODE	  	  
      !CHARACTER (LEN=500) :: IC_REF				
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbdiff
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbdiff2
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: EQMfname_PRN, VEQfname_PRN				
      CHARACTER (LEN=100) :: mesg

      INTEGER (KIND = prec_int4) :: J
      DOUBLE PRECISION MJDD0, MJDD, MJDref  
! ----------------------------------------------------------------------
	  REAL (KIND = prec_d) :: mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC
	  REAL (KIND = prec_d) :: dt_TT_TAI, dt_TAI_UTC, dt_TAI_GPS
      REAL (KIND = prec_d) :: t_sec     
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Delete Planetary ephemeris written file DE.430
fname_write = 'DE.430'
OPEN  (UNIT=7, FILE=fname_write)
CLOSE (UNIT=7, STATUS="DELETE")
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Satellite Orbits :: common configuration :: Forces model
! ----------------------------------------------------------------------
! Data reading: Gravitational Effects
! ----------------------------------------------------------------------
! General orbit parameterization											
!Call prm_main (EQMfname)
! Earth Gravity Field model
!CALL prm_gravity (EQMfname)												
! Planetary/Lunar ephemeris DE data 
!CALL prm_planets (EQMfname)												
! Ocean Tides model
!CALL prm_ocean (EQMfname)												
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellites Orbits :: PRN numbers 
! ----------------------------------------------------------------------
IF (IC_MODE_cfg == 1) THEN
param_id = 'pseudobs_filename'
CALL readparam (EQMfname, param_id, param_value)
ORBpseudobs_fname = param_value
!print *,"param_value", param_value 
!print *,"ORBpseudobs_fname", ORBpseudobs_fname 

! Read the sp3 header file :: Nsat
CALL sp3_PRN (ORBpseudobs_fname, PRNmatrix, Iyear, Imonth, Iday, Sec_00)
Nsat = size(PRNmatrix, DIM = 1)

ELSE IF (IC_MODE_cfg == 2) THEN

! Initial Conditions file 
CALL orbitIC (IC_filename_cfg, IC_matrix_glb, PRNmatrix)

Nsat = size(PRNmatrix, DIM = 1)

!print *,"IC_matrix_glb(1,1)", IC_matrix_glb(1,1)

mjd0   = IC_matrix_glb(1,1)
Sec_00 = IC_matrix_glb(1,2)

jd0 = 2400000.5D0
mjd = mjd0 + Sec_00 / 86400.0D0
CALL iau_JD2CAL ( jd0, mjd, Iyear, Imonth, Iday, FD, J_flag )

END IF
! ----------------------------------------------------------------------
print *,"Satellites number: ", Nsat, "IC Eopch: ", Iyear, Imonth, Iday, Sec_00
print *," "
! ----------------------------------------------------------------------
! Rewrite :: Initial Epoch
! ----------------------------------------------------------------------
! EQM & VEQ files
fname_id = '1'

param_id = 'Year'
write (param_value, *) Iyear
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'Month'
write (param_value, *) Imonth
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'Day'
write (param_value, *) Iday
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'Seconds'
write (param_value, *) Sec_00
! SCM 20190604 allow second > 10 to be written  - write (param_value, FMT='(F19.17)') Sec_00
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! Compute day of year
CALL iau_CAL2JD ( Iyear, Imonth, Iday, MJDD0, MJDD, J )
CALL iau_CAL2JD ( Iyear, 1, 1, MJDD0, MJDref, J )
DOY = IDNINT(MJDD-MJDref) + 1
YR = Iyear
!PRINT*,'Day Of Year =', Iyear,DOY

! Last modified: 19/02/2020 Thomas Papanikolaou: Correct MJD epoch for the computation of the time-variable gravity coefficients					
! ----------------------------------------------------------------------
! Satellite Orbits :: common configuration :: Forces model
! ----------------------------------------------------------------------
! Data reading: Gravitational Effects
! ----------------------------------------------------------------------
! General orbit parameterization											
Call prm_main (EQMfname)
! Earth Gravity Field model
CALL prm_gravity (EQMfname)												
! Planetary/Lunar ephemeris DE data 
CALL prm_planets (EQMfname)												
! Ocean Tides model
CALL prm_ocean (EQMfname)												
! ----------------------------------------------------------------------
!PRINT*,'NUMBER OF FORCE PARAMETERS =', NPARAM_glb
write(mesg, *) "A priori SRP model = ",SRP_MOD_arp
call report ('STATUS',pgrm_name,'',trim(mesg), '', 0)
if      (ECOM_param_glb /= 0 .and.  EMP_param_glb == 0) then
  write(mesg, *) "SRP model = ",ECOM_param_glb
else if (ECOM_param_glb == 0 .and.  EMP_param_glb /= 0) then
  write(mesg, *) "Empirical force model = ",EMP_param_glb
else if (ECOM_param_glb == 0 .and. EMP_param_glb == 0) then
  write(mesg, *) "WARNING: No SRP or Empirical force model estimated",EMP_param_glb,ECOM_param_glb
else if (ECOM_param_glb /= 0 .and. EMP_param_glb /= 0) then
  write(mesg, *) "Estimating EMP force model and ECOM SRP parameters estimated together: ",EMP_param_glb,ECOM_param_glb
else
  write(mesg, *) "Estimating EMP force model and ECOM SRP parameters not yet not supported: ",EMP_param_glb,ECOM_param_glb
  call report ('FATAL',pgrm_name,'pod_gnss',mesg, '/src/m_pod_gnss.f03', 1)
endif
call report ('STATUS',pgrm_name,'',mesg, ' ', 0)
print*,' '
! ----------------------------------------------------------------------
! Precise Orbit Determination :: Multi-GNSS multi-satellites POD loop
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
Do isat = 1 , Nsat

! ----------------------------------------------------------------------
! Modify/Rewrite the Configuration files
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite :: PRN
! ----------------------------------------------------------------------
PRN_isat = PRNmatrix(isat)
!print *,"Satellite: ", PRNmatrix(isat) ! isat
! Read Satellite information from SINEX file
! ----------------------------------------------------------------------
CALL read_satsnx (satsinex_filename_cfg, Iyear, DOY, Sec_00, PRN_isat)

write(*,10) trim(PRN_isat),SVNID,trim(BLKTYP),BLKID,POWER,MASS
10 format(' PRN: ',a,', SVN: ',i03,', BLK TYP: ',a,', BLKID: ',i3,', TX PWR: ',i3,', MASS: ',f8.3)
! ----------------------------------------------------------------------
! Copy Initial Configuration files 
write (fname_id, FMT='(A1,A3)') '_', PRN_isat
CALL write_prmfile2 (EQMfname, fname_id, EQMfname_PRN)
CALL write_prmfile2 (VEQfname, fname_id, VEQfname_PRN)
! ----------------------------------------------------------------------

write (fname_id, *) '_imd' !isat
param_id = 'Satellite_PRN'
param_value = PRN_isat
Call write_prmfile (EQMfname_PRN, fname_id, param_id, param_value)
Call write_prmfile (VEQfname_PRN, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite :: Initial State vector
! ----------------------------------------------------------------------
IF (IC_MODE_cfg == 1) THEN
! Interpolated Orbit: Read sp3 orbit data and apply Lagrange interpolation
Call prm_main     (EQMfname_PRN)
CALL prm_pseudobs (EQMfname_PRN)
Zo = pseudobs_ITRF(1,3:8)
print *,"IC: ", pseudobs_ITRF(1,:)

ELSE IF (IC_MODE_cfg == 2) THEN
! Initial Conditions file option
Zo = IC_matrix_glb (isat,3:8)
! Initial Conditions matrix per satellite
sz1 = size(IC_matrix_glb, DIM = 1)
sz2 = size(IC_matrix_glb, DIM = 2)
ALLOCATE (IC_sat_glb(sz2), STAT = AllocateStatus)
IC_sat_glb = IC_matrix_glb (isat,1:sz2)
print *,"Zo", IC_matrix_glb (isat,1:)
ELSE
print *,"Zo", Zo
END IF

! Write Initial Conditions (state vector only) in the configuration files
write (fname_id, *) '_imd' !isat
param_id = 'state_vector'
write (param_value, *) Zo
Call write_prmfile (EQMfname_PRN, fname_id, param_id, param_value)
Call write_prmfile (VEQfname_PRN, fname_id, param_id, param_value)
! Write Initial Conditions Reference System in the configuration files
IF (IC_MODE_cfg == 2) THEN
	write (fname_id, *) '_imd' !isat
	param_id = 'Reference_frame'
	write (param_value, *) IC_REF_cfg
	Call write_prmfile (EQMfname_PRN, fname_id, param_id, param_value)
	Call write_prmfile (VEQfname_PRN, fname_id, param_id, param_value)
END IF
! ----------------------------------------------------------------------

! End of update/rewrite configuration files
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Precise Orbit Determination :: main subroutine
!CAll orbitmain (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
CALL orbitmain (EQMfname_PRN, VEQfname_PRN, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms, &
		dorb_icrf, dorb_RTN, dorb_Kepler, dorb_itrf, orbdiff) 
! ----------------------------------------------------------------------
print *," "
print *," "


! ----------------------------------------------------------------------
! Allocation of the orbits & partial derivatives matrix
! ----------------------------------------------------------------------
if (isat == 1) then
sz1 = size(orb_icrf, DIM = 1)
sz2 = size(orb_icrf, DIM = 2)
Nepochs = sz1
N2_orb = sz2
!print *,"Nepochs: ", Nepochs

sz1 = size(veqSmatrix, DIM = 1)
sz2 = size(veqSmatrix, DIM = 2)
N2_veqSmatrix = sz2
!print *,"veqSmatrix n1: ", sz1

sz1 = size(veqPmatrix, DIM = 1)
sz2 = size(veqPmatrix, DIM = 2)
N2_veqPmatrix = sz2
!print *,"veqPmatrix n1: ", sz1

N2sum = 2 + (N2_orb-2) + (N2_veqSmatrix-2) + (N2_veqPmatrix-2)
!N2ics = 2 + (N2_veqSmatrix-2)/6 + (N2_veqPmatrix-2)/6
N2ics = 2 + 6 + NPARAM_glb
!print*,'N2ics: ',N2ics

! ----------------------------------------------------------------------
! Orbits matrix in ICRF
ALLOCATE (orbits_partials_icrf(Nepochs, N2sum , Nsat), STAT = AllocateStatus)
orbits_partials_icrf = 0.0D0
! Orbits matrix in ITRF
ALLOCATE (orbits_partials_itrf(Nepochs, N2sum , Nsat), STAT = AllocateStatus)
orbits_partials_itrf = 0.0D0
ALLOCATE (orbits_ics_icrf(N2ics , Nsat), STAT = AllocateStatus)
orbits_ics_icrf = 0.0D0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit comparison/residuals matrices
sz1 = size(dorb_RTN, DIM = 1)
sz2 = size(dorb_RTN, DIM = 2)
Ncommon = sz1

ALLOCATE (orbit_resR(Ncommon, Nsat+2), STAT = AllocateStatus)
orbit_resR = 0.0D0
ALLOCATE (orbit_resT(Ncommon, Nsat+2), STAT = AllocateStatus)
orbit_resT = 0.0D0
ALLOCATE (orbit_resN(Ncommon, Nsat+2), STAT = AllocateStatus)
orbit_resN = 0.0D0

orbit_resR(:,1:2) = dorb_RTN(:,1:2)
orbit_resT(:,1:2) = dorb_RTN(:,1:2)
orbit_resN(:,1:2) = dorb_RTN(:,1:2)
! ----------------------------------------------------------------------

sz1 = size(orbdiff, DIM = 1)
sz2 = size(orbdiff, DIM = 2)

ALLOCATE (orbdiff2(Nsat, sz1, sz2), STAT = AllocateStatus)
orbdiff2=0.0d0

end if
! ----------------------------------------------------------------------

orbdiff2 (isat,:,:) = orbdiff(:,:)

! ----------------------------------------------------------------------
! Create Orbit IC's matrix :: Write estimates for Satellite(isat) SVEC_Zo_ESTIM and ECOM_accel_aposteriori
! ----------------------------------------------------------------------
!orbits_ics_icrf(1:2,isat) = orb_icrf(1,1:2)

! ----------------------------------------------------------------------
! IC :: Initial Epoch MJD and Sec since 00h 
! ----------------------------------------------------------------------
! Time Scale transformation for the initial epoch time argument
! Time scale change is applied in case that TIME_SCALE .NOT. Terrestrial Time
! TIME_SCALE: global variable in module mdl_param.f03
! ----------------------------------------------------------------------
If (TIME_SCALE == 'TT') Then
! MJD_to and SEC_to :: global variables in mdl_param.f03
orbits_ics_icrf(1,isat) = MJD_to
orbits_ics_icrf(2,isat) = SEC_to

!Else If (TIME_SCALE /= 'TT') Then
Else 
! MJD_to and SEC_to :: global variables in mdl_param.f03

! MJD in TT Time
mjd = MJD_to
! Seconds since 00h
t_sec = SEC_to
!print *, "t_sec ", t_sec

! Time scale: TT to GPS time
CALL time_TT (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)
Call time_TT_sec (mjd , dt_TT_TAI, dt_TAI_UTC, dt_TAI_GPS)
!print *,"dt_TT_TAI, dt_TAI_UTC, dt_TAI_GPS", dt_TT_TAI, dt_TAI_UTC, dt_TAI_GPS

! Test the TIME_SCALE global variable in mdl_param.f03	
If (TIME_SCALE == 'GPS') then
	mjd = mjd_GPS
	t_sec = t_sec - (dt_TT_TAI + dt_TAI_GPS)
Else if (TIME_SCALE == 'UTC') then
	mjd = mjd_UTC		
	t_sec = t_sec - (dt_TT_TAI + dt_TAI_UTC)
Else if (TIME_SCALE == 'TAI') then
	mjd = mjd_TAI
	t_sec = t_sec - (dt_TT_TAI)	
End If	
!print *, "TIME_SCALE ", TIME_SCALE
!print *, "(dt_TT_TAI + dt_TAI_GPS) ", (dt_TT_TAI + dt_TAI_GPS)
!print *, "t_sec ", t_sec

orbits_ics_icrf(1,isat) = mjd
orbits_ics_icrf(2,isat) = t_sec

End If
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! IC :: Initial State Vector
! ----------------------------------------------------------------------
!orbits_ics_icrf(3:8,isat) = SVEC_Zo_ESTIM
orbits_ics_icrf(3:8,isat) = SVEC_Zo
!print *,"SVEC_Zo_ESTIM ", SVEC_Zo_ESTIM
!print *,"SVEC_Zo       ", SVEC_Zo
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! IC :: Orbital parameters being estimated e.g. force emprical parameters or SRP parameters
! ----------------------------------------------------------------------
!orbits_ics_icrf(9:8+(NPARAM_glb),isat) = ECOM_accel_aposteriori !*1.0D9
IF (NPARAM_glb /= 0) THEN
! ECOM SRP parameters estimated
  if ( ECOM_param_glb > 0 .and. EMP_param_glb == 0) then
    orbits_ics_icrf(9:8+(ECOMNUM),isat) = ECOM_accel_glb
!    orbits_ics_icrf(9:8+(NPARAM_glb),isat) = ECOM_accel_glb  ! Correction : Write the ECOM parameters in all POD cases including orbit propagation without estimation (POD MODE cases: 3 & 4)
! Empirical SRP parameters estimated
  else if ( EMP_param_glb > 0 .and. ECOM_param_glb == 0) then
    orbits_ics_icrf(9:11,isat)  = Bias_accel_glb  ! Write the EMP bias parameters
    orbits_ics_icrf(12:13,isat) = CPR_CS_glb(1,:) ! Write the EMP Radial CS  parameters
    orbits_ics_icrf(14:15,isat) = CPR_CS_glb(2,:) ! Write the EMP Transverse CS parameters
    orbits_ics_icrf(16:17,isat) = CPR_CS_glb(3,:) ! Write the EMP Normal CS  parameters
  else if (ECOM_param_glb > 0 .and. EMP_param_glb > 0) then
    orbits_ics_icrf(9:11,isat)  = Bias_accel_glb  
    orbits_ics_icrf(12:13,isat) = CPR_CS_glb(1,:) 
    orbits_ics_icrf(14:15,isat) = CPR_CS_glb(2,:) 
    orbits_ics_icrf(16:17,isat) = CPR_CS_glb(3,:)
    orbits_ics_icrf(18:17+(ECOMNUM),isat) = ECOM_accel_glb

  end if
END IF
!write(*,fmt='(a3,1x,f14.4,f14.6,1x,15(d17.10,1x))') PRN_isat,orbits_ics_icrf(:,isat)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit & Partial Derivatives matrix :: Write Orbit/Partials results from Satellite(isat) 
! ----------------------------------------------------------------------
!orbits_partials_icrf(:, 1:N2_orb , isat) = orb_icrf(:,:)
!orbits_partials_icrf(:, N2_orb + 1 : N2_orb + N2_veqSmatrix-2 , isat) = veqSmatrix(: , 3:N2_veqSmatrix)
!orbits_partials_icrf(:, N2_orb+N2_veqSmatrix-2 + 1 : N2_orb+N2_veqSmatrix-2 + N2_veqPmatrix-2 , isat)  &
!          = veqPmatrix(: , 3:N2_veqPmatrix)
		  
Do iepoch = 1 , Nepochs
Do iparam = 1 , N2_orb
orbits_partials_icrf(iepoch, iparam , isat) = orb_icrf(iepoch,iparam)
orbits_partials_itrf(iepoch, iparam , isat) = orb_itrf(iepoch,iparam)
End Do
End Do

Do iepoch = 1 , Nepochs
Do iparam = 1 , N2_veqSmatrix-2
orbits_partials_icrf(iepoch, N2_orb+iparam , isat) = veqSmatrix(iepoch , iparam+2)
orbits_partials_itrf(iepoch, N2_orb+iparam , isat) = veqSmatrix(iepoch , iparam+2)
End Do
End Do

Do iepoch = 1 , Nepochs
Do iparam = 1 , N2_veqPmatrix-2
orbits_partials_icrf(iepoch, N2_orb+N2_veqSmatrix-2+iparam , isat) = veqPmatrix(iepoch , iparam+2)
orbits_partials_itrf(iepoch, N2_orb+N2_veqSmatrix-2+iparam , isat) = veqPmatrix(iepoch , iparam+2)
End Do
End Do
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit residuals matrices per orbital component
! ----------------------------------------------------------------------
Do iepoch = 1 , Ncommon
orbit_resR(:,isat+2) = dorb_RTN(:,3)
orbit_resT(:,isat+2) = dorb_RTN(:,4)
orbit_resN(:,isat+2) = dorb_RTN(:,5)
End Do
! ----------------------------------------------------------------------
End Do


End SUBROUTINE

End MODULE

