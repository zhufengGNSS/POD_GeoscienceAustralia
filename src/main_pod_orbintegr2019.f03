      program main_pod


! ----------------------------------------------------------------------
! Program:	main_pod.f03
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
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      CHARACTER (LEN=100) :: EQMfname, VEQfname, PODfname, ORBMODfname				
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
      REAL (KIND = prec_d), DIMENSION(:,:,:), ALLOCATABLE :: orbit_veq  
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
fname_sp3 = 'igs16403.sp3'
! ----------------------------------------------------------------------
print *,"fname_sp3            ", fname_sp3
print *,"                     "


! ----------------------------------------------------------------------
! Configuration files of Orbit parameterization:
! ----------------------------------------------------------------------
EQMfname = 'EQM.in'
VEQfname = 'VEQ.in'
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

! ----------------------------------------------------------------------
! Temporary :: End of input POD configuration
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------  ! 999999999999999999999999999



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
!print *,"PRNmatrix  ", PRNmatrix
!print *,"PRNmatrix1 ", PRNmatrix(1)
!print *,"PRNmatrix2 ", PRNmatrix(2)


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


! ----------------------------------------------------------------------
! Precise Orbit Determination :: main subroutine
CAll orbitmain (EQMfname, VEQfname, orb_icrf, orb_itrf, veqSmatrix, veqPmatrix, Vres, Vrms)
! ----------------------------------------------------------------------

print *," "
print *," "


! ----------------------------------------------------------------------
! Allocation of the orbit & partial derivatives matrix
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

ALLOCATE (orbit_veq(Nepochs, N2sum , Nsat), STAT = AllocateStatus)
orbit_veq = 0.0D0


ALLOCATE (orbit_resR(Nepochs, Nsat+2), STAT = AllocateStatus)
orbit_resR = 0.0D0
ALLOCATE (orbit_resT(Nepochs, Nsat+2), STAT = AllocateStatus)
orbit_resT = 0.0D0
ALLOCATE (orbit_resN(Nepochs, Nsat+2), STAT = AllocateStatus)
orbit_resN = 0.0D0

orbit_resR(:,1:2) = Vres(:,1:2)
orbit_resT(:,1:2) = Vres(:,1:2)
orbit_resN(:,1:2) = Vres(:,1:2)

!print *,"isat", isat
!print *,"orbit_veq", orbit_veq
!print *,"Nepochs", Nepochs
!print *,"N2_orb", N2_orb
!print *,"N2_veqSmatrix", N2_veqSmatrix
!print *,"N2_veqPmatrix", N2_veqPmatrix
!print *,"N2sum", N2sum
end if
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Orbit & Partial Derivatives matrix :: Write Orbit/Partials results from Satellite(isat) 
! ----------------------------------------------------------------------
!orbit_veq(:, 1:N2_orb , isat) = orb_icrf(:,:)
!orbit_veq(:, N2_orb + 1 : N2_orb + N2_veqSmatrix-2 , isat) = veqSmatrix(: , 3:N2_veqSmatrix)
!orbit_veq(:, N2_orb+N2_veqSmatrix-2 + 1 : N2_orb+N2_veqSmatrix-2 + N2_veqPmatrix-2 , isat)  &
!          = veqPmatrix(: , 3:N2_veqPmatrix)
		  
Do iepoch = 1 , Nepochs
Do iparam = 1 , N2_orb
orbit_veq(iepoch, iparam , isat) = orb_icrf(iepoch,iparam)
End Do
End Do

Do iepoch = 1 , Nepochs
Do iparam = 1 , N2_veqSmatrix-2
orbit_veq(iepoch, N2_orb+iparam , isat) = veqSmatrix(iepoch , iparam+2)
End Do
End Do

Do iepoch = 1 , Nepochs
Do iparam = 1 , N2_veqPmatrix-2
orbit_veq(iepoch, N2_orb+N2_veqSmatrix-2+iparam , isat) = veqPmatrix(iepoch , iparam+2)
End Do
End Do
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit residuals matrices per orbital component
! ----------------------------------------------------------------------
Do iepoch = 1 , Nepochs

orbit_resR(:,isat+2) = Vres(:,3)
orbit_resT(:,isat+2) = Vres(:,4)
orbit_resN(:,isat+2) = Vres(:,5)

End Do
! ----------------------------------------------------------------------

End Do


! ----------------------------------------------------------------------
! Write satellite orbits and partial derivatives to one .orb output file (internal format)
! ----------------------------------------------------------------------
orbits_fname = 'orbits_partials.orb'
CALL writeorbit_multi (orbit_veq, PRNmatrix, orbits_fname)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Write satellite orbits to sp3 format
! ----------------------------------------------------------------------
! CALL m_writesp3 (orbit_veq, PRNmatrix, sp3filename, sat_prn, sat_vel)
ORB2sp3_fname = 'GA.sp3'
!write (ORB2sp3_fname, *) 'GA', '' ,'.sp3'
CALL write_orb2sp3 (orbit_veq, PRNmatrix, ORB2sp3_fname, sat_vel)
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


!isat = 1
!print *,"orbit_veq      ", orbit_veq
!print *,"orbit_veq 1:1  ", orbit_veq(1,:,1)
!print *,"orbit_veq 111  ", orbit_veq(1,1,isat), isat
!print *,"orbit_veq 121  ", orbit_veq(1,2,isat)
!print *,"orbit_veq 131  ", orbit_veq(1,3,isat)
!print *,"orbit_veq 211  ", orbit_veq(2,1,isat), isat
!print *,"orbit_veq ", orbit_veq(1,1:2,isat)
!print *,"orbit_veq ", orbit_veq(1,3:8,isat)
!print *,""

!Do i = 1 , N2sum
!print *,"orbit_veq( ,i,)", orbit_veq(1,i,1)
!print *,"orbit_veq( ,i,)", orbit_veq(2,i,1)
!End Do
!print *,""

!print *,"orb_icrf  ", orb_icrf(1,:)
!print *,""
!print *,"veqSmatrix", veqSmatrix(1,:)
!print *,""
!print *,"veqPmatrix", veqPmatrix(1,:)

!sz1 = size(orbit_veq, DIM = 1)
!print *,"sz1", sz1
!sz1 = size(orbit_veq, DIM = 2)
!print *,"sz1", sz1
!sz1 = size(orbit_veq, DIM = 3)
!print *,"sz1",sz1


CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


End Program

