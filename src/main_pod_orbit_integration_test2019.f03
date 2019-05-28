      program main_pod_orbit_integration_test2019


! ----------------------------------------------------------------------
! Program:	main_pod_orbit_integration_test2019.f03
! ----------------------------------------------------------------------
! Purpose:
!  Precise Orbit Determination (POD) of GNSS satellites 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, CRC-SI
! Created:	13 September 2017
! ----------------------------------------------------------------------
! POD version major modifications highlights: 
! Last modified  
! - Dr. Thomas Papanikolaou, 3 May 2018
! 	Preliminary version of GNSS dynamic orbit determination	
! - Dr. Thomas Papanikolaou, 25 June 2018
! 	Version with minor revisions
! - Dr. Thomas Papanikolaou, 30 November 2018
! 	Precise Orbit Determination (POD) version: Estimation of empirical forces parameters (bias, cycle-per-rev) that lead to mm-cm level orbital accuracy w.r.t. IGS precise orbits
! - Dr. Thomas Papanikolaou, 30 January 2019
! 	POD version upgrade: Ocean tides effect revision that has significant impact on longer orbit arcs e.g. 3 days 
! - Dr. Thomas Papanikolaou, 29 March 2019
! 	POD version upgrade to a multi-GNSS multi-satellite POD version 
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
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
	  USE m_orbext2
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: EQMfname, VEQfname, PODfname, ORBMODfname				
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orb_icrf, orb_itrf  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqSmatrix, veqPmatrix
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Vres  
      REAL (KIND = prec_d), DIMENSION(3) :: Vrms 	    
	  REAL (KIND = prec_d), DIMENSION(5,6) :: stat_XYZ_extC, stat_RTN_extC, stat_Kepler_extC, stat_XYZ_extT
! ----------------------------------------------------------------------
      CHARACTER (LEN=2) :: GNSS_id
	  INTEGER (KIND = prec_int2) :: ORB_mode
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Nsat, isat
	  INTEGER (KIND = prec_int8) :: iepoch, iparam
	  INTEGER (KIND = prec_int8) :: i
	  INTEGER (KIND = prec_int8) :: sz1, sz2, Nepochs, N2_orb, N2_veqSmatrix, N2_veqPmatrix, N2sum  
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_icrf  
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_itrf  
	  CHARACTER (LEN=3), ALLOCATABLE :: PRNmatrix(:)
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus  
	  CHARACTER (LEN=3) :: PRN_isat
	  INTEGER :: ios
      CHARACTER (LEN=100) :: orbits_fname				
      CHARACTER (LEN=100) :: fname_write				
      CHARACTER (LEN=100) :: filename				
! ----------------------------------------------------------------------
      CHARACTER (LEN=300) :: fname_sp3				
	  INTEGER :: year, month, day
	  INTEGER :: Iyear, Imonth, Iday
      REAL (KIND = prec_d) :: Sec_00 	    
! ----------------------------------------------------------------------
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
      REAL (KIND = prec_d) :: Zo(6) 
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: ORB2sp3_fname				
      INTEGER (KIND = prec_int2) :: sat_vel	  	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resR  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resT  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resN  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: mjd
      INTEGER (KIND = prec_int8) :: GPS_week, GPSweek_mod1024
      REAL (KIND = prec_d) :: GPS_wsec, GPS_day
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: force_model, institute_sp3	  	  
      CHARACTER (LEN=100) :: EQMfname_G01, EQMfname_G02
      CHARACTER (LEN=100) :: fname_sp3_eci, fname_sp3_ecf
! ----------------------------------------------------------------------
	  INTEGER (KIND = prec_int8) :: Ncommon  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_icrf, dorb_itrf 
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: dorb_RTN, dorb_Kepler
! ----------------------------------------------------------------------
      CHARACTER (LEN=11) :: GFM_fname
      CHARACTER (LEN=1) :: pole_tide
      !INTEGER (KIND = prec_int2) :: pole_tide	  	  




! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit modelling cases :: 
! ----------------------------------------------------------------------
! 1. not_nst
! 2. not
! 3. wot
force_model = 1
! ----------------------------------------------------------------------
! Orbit comparison :: Institute sp3 files
! ----------------------------------------------------------------------
! 1. ESOC
! 2. JPL
! 3. MIT
! 4. TUG
institute_sp3 = 3
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Gravity Field model
GFM_fname = 'goco05s.gfc'
GFM_fname = 'egm2008.gfc'
! ----------------------------------------------------------------------
! Pole Tide (Solid Earth Tide & Ocean Pole tide)
! 0. not included 
! 1. included
pole_tide = '0'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------




! CPU Time
CALL cpu_time (CPU_t0)


! ----------------------------------------------------------------------  ! 999999999999999999999999999
! ----------------------------------------------------------------------
! Temporary :: (To be upgraded to one major input configuration file)
! ----------------------------------------------------------------------
! Start of input POD configuration
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Input a-priori orbit in sp3 format (to be used as pseudo-observations)
!fname_sp3 = 'igs16403.sp3'
!fname_sp3 = 'igs19423.sp3'
!fname_sp3 = 'igs19424.sp3'
! ----------------------------------------------------------------------
!print *,"fname_sp3            ", fname_sp3
!print *,"                     "


! ----------------------------------------------------------------------
! Configuration files of Orbit parameterisation:
! ----------------------------------------------------------------------
!EQMfname = 'EQM.in'
!VEQfname = 'VEQ.in'

! ----------------------------------------------------------------------
IF (force_model == 1) THEN

EQMfname_G01 = 'IGS_test_G01_EQM_not_nst.in'
EQMfname_G02 = 'IGS_test_G02_EQM_not_nst.in'
! ----------------------------------------------------------------------
IF (institute_sp3 == 1) THEN
fname_sp3_eci = 'esoc2_eci_not_nst.sp3'
fname_sp3_ecf = 'esoc2_ecf_not_nst.sp3'
ELSE IF (institute_sp3 == 2) THEN
fname_sp3_eci = 'jpl_eci_not_nst.sp3'
fname_sp3_ecf = 'jpl_ecf_not_nst.sp3'
ELSE IF (institute_sp3 == 3) THEN
fname_sp3_eci = 'mit_eci_not_nst.sp3'
fname_sp3_ecf = 'mit_ecf_not_nst.sp3'
ELSE IF (institute_sp3 == 4) THEN
fname_sp3_eci = 'tug_eci_not_nst.sp3'
fname_sp3_ecf = 'tug_ecf_not_nst.sp3'
END IF
! ----------------------------------------------------------------------

ELSE IF (force_model == 2) THEN

EQMfname_G01 = 'IGS_test_G01_EQM_not.in'
EQMfname_G02 = 'IGS_test_G02_EQM_not.in'
! ----------------------------------------------------------------------
IF (institute_sp3 == 1) THEN
fname_sp3_eci = 'esoc2_eci_not.sp3'
fname_sp3_ecf = 'esoc2_ecf_not.sp3'
ELSE IF (institute_sp3 == 2) THEN
fname_sp3_eci = 'jpl_eci_not.sp3'
fname_sp3_ecf = 'jpl_ecf_not.sp3'
ELSE IF (institute_sp3 == 3) THEN
fname_sp3_eci = 'mit_eci_not.sp3'
fname_sp3_ecf = 'mit_ecf_not.sp3'
ELSE IF (institute_sp3 == 4) THEN
fname_sp3_eci = 'tug_eci_not.sp3'
fname_sp3_ecf = 'tug_ecf_not.sp3'
END IF
! ----------------------------------------------------------------------

ELSE IF (force_model == 3) THEN

EQMfname_G01 = 'IGS_test_G01_EQM_wot.in'
EQMfname_G02 = 'IGS_test_G02_EQM_wot.in'
! ----------------------------------------------------------------------
IF (institute_sp3 == 1) THEN
fname_sp3_eci = 'esoc2_eci_wot.sp3'
fname_sp3_ecf = 'esoc2_ecf_wot.sp3'
ELSE IF (institute_sp3 == 2) THEN
fname_sp3_eci = 'jpl_eci_wot.sp3'
fname_sp3_ecf = 'jpl_ecf_wot.sp3'
ELSE IF (institute_sp3 == 3) THEN
fname_sp3_eci = 'mit_eci_wot.sp3'
fname_sp3_ecf = 'mit_ecf_wot.sp3'
ELSE IF (institute_sp3 == 4) THEN
fname_sp3_eci = 'tug_eci_wot.sp3'
fname_sp3_ecf = 'tug_ecf_wot.sp3'
END IF
! ----------------------------------------------------------------------

END IF
! ----------------------------------------------------------------------
EQMfname = EQMfname_G01
!VEQfname = 'VEQ.in'
VEQfname = EQMfname_G01

! Read the PRNs for PRN_matrix
fname_sp3 = fname_sp3_eci

print *,"force_model   ", force_model
print *,"EQMfname      ", EQMfname
print *,"VEQfname      ", VEQfname
print *,"fname_sp3_eci ", fname_sp3_eci
print *,"fname_sp3_ecf ", fname_sp3_ecf
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Temporary:: manual configuration of 3 global parameters
! ----------------------------------------------------------------------
! "GPS satellite Block" and "Beidou Orbit Type"
! ----------------------------------------------------------------------
! GPS case: Satellite Block ID:	1=I, 2=II, 3=IIA, IIR=(4, 5), IIF=6
SATblock_glb = 4
! ----------------------------------------------------------------------
! Beidou case:
! 1. BDSorbtype = 'IGSO'
! 2. BDSorbtype = 'MEO'
BDSorbtype_glb = 'MEO'
! ----------------------------------------------------------------------
! Empirical Forces reference frame:
! 1. Orbital Frame
! 2. Body-fixed frame
Frame_EmpiricalForces_glb = 1
! ----------------------------------------------------------------------
print *,"              "
print *,"Frame_EmpiricalForces_glb ", Frame_EmpiricalForces_glb
print *,"SATblock_glb              ", SATblock_glb
print *,"BDSorbtype_glb            ", BDSorbtype_glb
print *,"              "

! ----------------------------------------------------------------------
! Write option for Satellite Velocity vector in computed orbit to output sp3 format 
! ----------------------------------------------------------------------
! 0. sat_vel = 0 :: Do not write Velocity vector to sp3 orbit
! 1. sat_vel > 0 :: Write Velocity vector to sp3 orbit
sat_vel = 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Temporary :: End of input POD configuration
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------  ! 999999999999999999999999999


! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set a-priori orbit (pseudo-observations; orbit comparison)
! ----------------------------------------------------------------------
!write (fname_id, *) isat
fname_id = '1'

param_id = 'pseudobs_filename'
param_value = fname_sp3_ecf
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'orbit_filename'
param_value = fname_sp3_eci
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Rewrite Configuration files :: Modify the gravity field model
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'gravity_model_filename'
param_value = GFM_fname
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Rewrite Configuration files :: Pole Tide effect
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'solid_earth_pole_tide'
param_value = pole_tide
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

fname_id = '1'
param_id = 'ocean_pole_tide'
param_value = pole_tide
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
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
Call prm_main (EQMfname)
! Earth Gravity Field model
CALL prm_gravity (EQMfname)												
! Planetary/Lunar ephemeris DE data 
CALL prm_planets (EQMfname)												
! Ocean Tides model
CALL prm_ocean (EQMfname)												
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellites Orbits :: PRN numbers 
! ----------------------------------------------------------------------
! Read the sp3 header file :: Nsat
CALL sp3_PRN (fname_sp3, PRNmatrix, Iyear, Imonth, Iday, Sec_00)
Nsat = size(PRNmatrix, DIM = 1)
! ----------------------------------------------------------------------
print *,"Satellites Orbits number: ", Nsat
print *," "





! ---------------------------------------------------------------------- ! Deactivate
if (1<0) then
! ----------------------------------------------------------------------
! Rewrite :: Initial Epoch
! ----------------------------------------------------------------------
write (fname_id, *) Nsat
! ----------------------------------------------------------------------
! EQM file
ORBMODfname = EQMfname

param_id = 'Year'
write (param_value, *) Iyear
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)

param_id = 'Month'
write (param_value, *) Imonth
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)

param_id = 'Day'
write (param_value, *) Iday
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)

param_id = 'Seconds'
write (param_value, FMT='(F19.17)') Sec_00
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------
! VEQ file
ORBMODfname = VEQfname

param_id = 'Year'
write (param_value, *) Iyear
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)

param_id = 'Month'
write (param_value, *) Imonth
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)

param_id = 'Day'
write (param_value, *) Iday
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)

param_id = 'Seconds'
write (param_value, FMT='(F19.17)') Sec_00
Call write_prmfile (ORBMODfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------
!print *,"param_id ", param_id
!print *,"param_value ", param_value
end if
! ---------------------------------------------------------------------- ! Deactivate





! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Precise Orbit Determination :: Multi-GNSS multi-satellites POD loop
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
Do isat = 1 , 1 ! Nsat



! ---------------------------------------------------------------------- ! Deactivate
if (1<0) then 
! ----------------------------------------------------------------------
! Modify/Rewrite the Configuration files
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite :: PRN
! ----------------------------------------------------------------------
PRN_isat = PRNmatrix(isat)
print *,"Satellite: ", PRNmatrix(isat) ! isat
!print *,"PRN       ",PRN

write (fname_id, *) isat
param_id = 'Satellite_PRN'
!write (param_value, *) PRN_isat
param_value = PRN_isat
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Rewrite :: Initial State vector
! ----------------------------------------------------------------------
! Interpolated Orbit: Read sp3 orbit data and apply Lagrange interpolation
Call prm_main     (EQMfname)
CALL prm_pseudobs (EQMfname)
Zo = pseudobs_ITRF(1,3:8)
!print *,"Zo", Zo

write (fname_id, *) isat
param_id = 'state_vector'
write (param_value, *) Zo
!param_value = PRN_isat
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

!INTEGER (KIND = prec_int8) :: interv, NPint
!REAL (KIND = prec_q), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: orbint
!CALL interp_orb (fname_sp3, PRN_isat, interpstep, NPint, orbsp3)
! Initial State Vector
! Zo = orbsp3 (1,3:8) 
! ----------------------------------------------------------------------

! End of update/rewrite configuration files
! ----------------------------------------------------------------------
end if
! ---------------------------------------------------------------------- ! Deactivate



! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
print *,"isat", isat

IF (isat == 1) THEN

EQMfname = EQMfname_G01
print *,"EQMfname   ", EQMfname
print *,"Satellite: ", PRNmatrix(isat) ! isat
print *," "

ELSE IF (isat == 2) THEN

EQMfname = EQMfname_G02
print *,"EQMfname   ", EQMfname
print *,"Satellite: ", PRNmatrix(isat) ! isat
print *," "

END IF


! ----------------------------------------------------------------------
! Rewrite Configuration files :: Set sp3 file name for orbit comparison)
! ----------------------------------------------------------------------
!write (fname_id, *) isat
fname_id = '1'

param_id = 'pseudobs_filename'
param_value = fname_sp3_ecf
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

param_id = 'orbit_filename'
param_value = fname_sp3_eci
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Rewrite Configuration files :: Modify the gravity field model
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'gravity_model_filename'
param_value = GFM_fname
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Rewrite Configuration files :: Pole Tide effect
! ----------------------------------------------------------------------
fname_id = '1'
param_id = 'solid_earth_pole_tide'
param_value = pole_tide
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)

fname_id = '1'
param_id = 'ocean_pole_tide'
param_value = pole_tide
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Cases of JPL, MIT :: sp3 orbits without Velocity vector
! ----------------------------------------------------------------------
IF (institute_sp3 == 2 .or. institute_sp3 == 3) THEN
fname_id = '1'
param_id = 'orbit_external_opt'
write (param_value, *) 2
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
END IF
! ----------------------------------------------------------------------

 
! ----------------------------------------------------------------------
! Precise Orbit Determination :: main subroutine
!CAll orbitmain (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
CALL orbitmain (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms, &
				dorb_icrf, dorb_RTN, dorb_Kepler, dorb_itrf) 
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

sz1 = size(veqSmatrix, DIM = 1)
sz2 = size(veqSmatrix, DIM = 2)
N2_veqSmatrix = sz2

sz1 = size(veqPmatrix, DIM = 1)
sz2 = size(veqPmatrix, DIM = 2)
N2_veqPmatrix = sz2

N2sum = 2 + (N2_orb-2) + (N2_veqSmatrix-2) + (N2_veqPmatrix-2)

! ----------------------------------------------------------------------
! Orbits matrix in ICRF
ALLOCATE (orbits_partials_icrf(Nepochs, N2sum , Nsat), STAT = AllocateStatus)
orbits_partials_icrf = 0.0D0
! Orbits matrix in ITRF
ALLOCATE (orbits_partials_itrf(Nepochs, N2sum , Nsat), STAT = AllocateStatus)
orbits_partials_itrf = 0.0D0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit Comparison Matrices
sz1 = size(dorb_RTN, DIM = 1)
sz2 = size(dorb_RTN, DIM = 2)
Ncommon = sz1

ALLOCATE (orbit_resR(Ncommon, Nsat+2), STAT = AllocateStatus)
orbit_resR = 0.0D0
ALLOCATE (orbit_resT(Ncommon, Nsat+2), STAT = AllocateStatus)
orbit_resT = 0.0D0
ALLOCATE (orbit_resN(Ncommon, Nsat+2), STAT = AllocateStatus)
orbit_resN = 0.0D0

!orbit_resR(:,1:2) = dorb_RTN(:,1:2)
!orbit_resT(:,1:2) = dorb_RTN(:,1:2)
!orbit_resN(:,1:2) = dorb_RTN(:,1:2)
! ----------------------------------------------------------------------
end if
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
! Orbit residuals/comparison matrices per orbital component
! ----------------------------------------------------------------------
Do iepoch = 1 , Ncommon
orbit_resR(iepoch,isat+2) = dorb_RTN(iepoch,3)
orbit_resT(iepoch,isat+2) = dorb_RTN(iepoch,4)
orbit_resN(iepoch,isat+2) = dorb_RTN(iepoch,5)
END DO
! ----------------------------------------------------------------------

End Do

orbit_resR(:,1:2) = dorb_RTN(:,1:2)
orbit_resT(:,1:2) = dorb_RTN(:,1:2)
orbit_resN(:,1:2) = dorb_RTN(:,1:2)
!Do iepoch = 1 , Ncommon
!orbit_resR(iepoch,1) = dorb_RTN(iepoch,1)
!orbit_resR(iepoch,2) = dorb_RTN(iepoch,2)
!orbit_resT(iepoch,1) = dorb_RTN(iepoch,1)
!orbit_resT(iepoch,2) = dorb_RTN(iepoch,2)
!orbit_resN(iepoch,1) = dorb_RTN(iepoch,1)
!orbit_resN(iepoch,2) = dorb_RTN(iepoch,2)
!END DO



! ----------------------------------------------------------------------
! Write satellite orbits and partial derivatives to one .orb output file (internal format)
! ----------------------------------------------------------------------
orbits_fname = 'orbits_partials_icrf.orb'
CALL writeorbit_multi (orbits_partials_icrf, PRNmatrix, orbits_fname)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Write satellite orbits to sp3 format
! ----------------------------------------------------------------------
! ICRF
ORB2sp3_fname = 'ga_eci.sp3'
CALL write_orb2sp3 (orbits_partials_icrf, PRNmatrix, ORB2sp3_fname, sat_vel)

! ITRF
ORB2sp3_fname = 'ga_ecf.sp3'
CALL write_orb2sp3 (orbits_partials_itrf, PRNmatrix, ORB2sp3_fname, sat_vel)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write Orbit residuals
! ----------------------------------------------------------------------
filename = "orbit_residuals_R.out"
Call writearray (orbit_resR, filename)
filename = "orbit_residuals_T.out"
Call writearray (orbit_resT, filename)
filename = "orbit_residuals_N.out"
Call writearray (orbit_resN, filename)
! ----------------------------------------------------------------------



CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0
PRINT *,"End Program"
CALL EXIT(0)

End 

