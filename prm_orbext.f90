SUBROUTINE prm_orbext (orbext)


! ----------------------------------------------------------------------
! SUBROUTINE: prm_orbext.f90
! ----------------------------------------------------------------------
! Purpose:
!  External Orbit settings
!  Precise orbit data (sp3) use or computation of Keplerian orbits is set here.
! 
!  The precise orbit data (sp3) are used for the following purposes:
!  1. External orbit comparison 
!  2. Pseudo-observations within the dynamic orbit estimation procedure 
! ----------------------------------------------------------------------
! Remark 1:
!  The option of computing a Keplerian orbit (instead of using orbit sp3 data) 
!  is supported and configured here
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Cooperative Research Centre for Spatial Information, Australia
! Created:	20 September 2017
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_legendre
      USE mdl_legendre1
      USE mdl_gfc
      USE mdl_planets
      USE mdl_tides	  
      USE mdl_write
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: data_opt
      CHARACTER (LEN=300) :: fname_orb
      CHARACTER (LEN=300) :: fname_orb_0, fname_orb_1, fname_orb_2
      CHARACTER (LEN=300) :: fname_orbint
      CHARACTER (LEN=300) :: fname_write
      INTEGER (KIND = prec_int8) :: NPint
      INTEGER (KIND = prec_int8) :: interpstep
      INTEGER (KIND = prec_int8) :: sz1, sz2 
      INTEGER (KIND = prec_int8) :: Ndays
      INTEGER (KIND = prec_int4) :: Zo_el
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! External Orbit
! ----------------------------------------------------------------------
! 1. RSO orbit data; GPS Position and Veloctiy vectors; Interval: 30 sec
! 2. Orbit based on Lagrange interpolation of sp3 data e.g. IGS final/rapid (15 min); MGEX (5 min)
! 3. Orbit: Position vector from input sp3 file; Velocity vector approximated through position differences (sp3) without applying interpolation 
! 4. Orbit arc (Lagrange interpolation) based on 3 boundary sp3 files 
! 5: Long orbit arc based on Keplerian orbit 
data_opt = orbext
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! GNSS orbit data (sp3) file name
fname_orb = 'igs18885.sp3'	! 18/03/2016  IIF
! Case 4:
!fname_orb_0 = 'wum19283.sp3'  
!fname_orb_1 = 'wum19284.sp3'  
!fname_orb_2 = 'wum19285.sp3'  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Interpolated Orbit:	Cases: 2,4
! ----------------------------------------------------------------------
! Interpolation interval
interpstep = 300 
!interpstep = integstep
! Number of data points used in Lagrange interpolation   
NPint = 12
! Output file name for writing interpolated orbit (Cases: 2, 4)
fname_orbint = 'orb_interp.out'	  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Keplerian orbit:		Case: 5 
! ----------------------------------------------------------------------
! Keplerian orbit arc length in number of days	  
Ndays = 352 ! GPS draconitic year 351.4d

! Initial state vector type:
! 1. Keplerian Elements (Inertial frame)
! 2. State Vector (Position & Velocity) cartesian coordinates (in ICRF)
Zo_el = 2	
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Read GNSS orbit data per epoch and form the orbit arrays (satellite state vector per epoch)
! ----------------------------------------------------------------------
 
! ----------------------------------------------------------------------
! Case 1: RSO data by GFZ; Data interval 30 sec
If (data_opt == 1) then 
      Call sp3 (fname_orb, PRN)
	  sz1 = size(orb1, DIM = 1)
	  sz2 = size(orb1, DIM = 2)
	  
      ALLOCATE ( orb2(sz1, sz2) , STAT = AllocateStatus)
      orb2 = orb1

! Write Orbit after interpolation 
      ALLOCATE (wrtArray( sz1, sz2), STAT = AllocateStatus)
      fname_write = 'orb_rso.out'
      wrtArray = orb2
      CALL write_array (fname_write)  
      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! Case 2: Orbit based on Lagrange interpolation of sp3 data e.g. IGS final/rapid (15 min); MGEX (5 min) IGS rapid orbits; interval 15 min
Else if (data_opt == 2) then
	  Call interp_orb(fname_orb, PRN, interpstep, NPint, fname_orbint)
	  sz1 = size(orb1, DIM = 1)
      sz2 = size(orb1, DIM = 2)
      ! Allocatable array for the orbit (position vector only) in ITRF
      ALLOCATE (orb3(sz1,sz2), STAT = AllocateStatus) 
	  orb3 = orb1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Case 3: Orbit: Position from input sp3 file; Velocity approximated through position differences without applying interpolation  
Else if (data_opt == 3) then
      fname_write = 'orb_sp3.out' 
	  CALL orb_sp3 (fname_orb, PRN, fname_write)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Case 4: Long orbit arc obtained from multiple sp3 files 
else if (data_opt == 4) then
	! Performing interpolation individually for each sp3 orbit file
	! Extrapolation is applied in order to avoid the gaps at the end of the daily file (sp3 last epoch is 23h 45min 00sec)
	CALL orb3sp3(fname_orb_0, fname_orb_1, fname_orb_2, PRN, interpstep, NPint, fname_orbint)	
! or	
	! Performing interpolation after forming the orbit array (position vector) based on the mutliple sp3 files
	!CALL orb3sp3_2(fname_orb_0, fname_orb_1, fname_orb_2, PRN, interpstep, NPint, fname_orbint)	
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Case 5: Long orbit arc during draconitic periods based on GNSS Keplerian orbits
else if (data_opt == 5) then	  
	  CALL keplerorb (MJDo, Sec0, Zo, Zo_el, Ndays, interpstep)

	  sz1 = size(orb1, DIM = 1)
	  sz2 = size(orb1, DIM = 2)
      ALLOCATE (wrtArray(sz1,sz2), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0) then 
		print *, " "
		print *, "*** Insufficient memory ***"
		print *, "*** Allocatable array: wrtArray for write_array.f90 : Writing Kepler orbit ***"
		print *, " "
	  End If

      fname_write = 'orb_kepler2.out'
      wrtArray = orb1
      CALL write_array (fname_write)
	  
      DEALLOCATE (wrtArray, STAT = DeAllocateStatus)
	  DEALLOCATE (orb1, STAT = DeAllocateStatus)	  
! ----------------------------------------------------------------------
	  
End if	  





END
