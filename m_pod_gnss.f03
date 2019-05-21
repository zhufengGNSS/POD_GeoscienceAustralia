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


SUBROUTINE pod_gnss (EQMfname, VEQfname, PRNmatrix, orbits_partials_icrf, orbits_partials_itrf, orbit_resR, orbit_resT, orbit_resN)

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
	  INTEGER (KIND = prec_int8) :: sz1, sz2, Nepochs, N2_orb, N2_veqSmatrix, N2_veqPmatrix, N2sum  
      !REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_icrf  
      !REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbits_partials_itrf  
	  !CHARACTER (LEN=3), ALLOCATABLE :: PRNmatrix(:)
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus  
	  CHARACTER (LEN=3) :: PRN_isat
	  INTEGER :: ios
      CHARACTER (LEN=100) :: orbits_fname				
      CHARACTER (LEN=100) :: fname_write				
      CHARACTER (LEN=100) :: filename				
! ----------------------------------------------------------------------
      CHARACTER (LEN=300) :: fname_sp3, ORBpseudobs_fname, ORBEXT_fname				
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
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resR  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resT  
      !REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: orbit_resN  
! ----------------------------------------------------------------------
      !REAL (KIND = prec_d) :: mjd
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
param_id = 'pseudobs_filename'
CALL readparam (EQMfname, param_id, param_value)
ORBpseudobs_fname = param_value
!print *,"param_value", param_value 
!print *,"ORBpseudobs_fname", ORBpseudobs_fname 

! Read the sp3 header file :: Nsat
CALL sp3_PRN (ORBpseudobs_fname, PRNmatrix, Iyear, Imonth, Iday, Sec_00)
Nsat = size(PRNmatrix, DIM = 1)
! ----------------------------------------------------------------------
print *,"Satellites number: ", Nsat
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
write (param_value, FMT='(F19.17)') Sec_00
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
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
print *,"Satellite: ", PRNmatrix(isat) ! isat

write (fname_id, *) isat
param_id = 'Satellite_PRN'
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
Call write_prmfile (EQMfname, fname_id, param_id, param_value)
Call write_prmfile (VEQfname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

! End of update/rewrite configuration files
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

