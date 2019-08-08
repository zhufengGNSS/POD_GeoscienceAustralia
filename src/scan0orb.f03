SUBROUTINE scan0orb 

! ----------------------------------------------------------------------
! SUBROUTINE: scan0orb.f03
! ----------------------------------------------------------------------
! Purpose:
!  To skip the bad orbit with zero value  in SP3 file
! ----------------------------------------------------------------------
! Input arguments:
! 
! - pseudobs_ICRF: 	Orbit array (Nx5) in inertial frame (ICRF) including the state vector per epoch
! 					Collumns desciption per epoch:
!               	- Modified Julian Day number (including the fraction of the day) 
!					- Seconds since 00h 
!					- Position vector (m)
! - pseudobs_ITRF:	Orbit array (Nx5) in terrestrial frame (ITRF)
! 					Collumns desciption per epoch: similar to the orbobs_ICRF
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng,  Geoscience Australia
!
! Created:	07-08-2019
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param

      IMPLICIT NONE

	  
 
! ----------------------------------------------------------------------
! Local Variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: sz1, sz2, sz3 ,sz4 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, read_i, k 
      INTEGER (KIND = prec_int8) :: ic, it, j
! ----------------------------------------------------------------------

sz1 = SIZE (pseudobs_ICRF, DIM =1)
sz2 = SIZE (pseudobs_ICRF, DIM =2)

sz3 = SIZE (pseudobs_ITRF, DIM =1)
sz4 = SIZE (pseudobs_ITRF, DIM =2)

!print*,'sz1, sz2 =', sz1, sz2, 'sz3, sz4 =',  sz3, sz4
pseudobs_ICRF2 = pseudobs_ICRF
pseudobs_ITRF2 = pseudobs_ITRF

pseudobs_ICRF2 = 0.d0
pseudobs_ITRF2 = 0.d0
! ----------------------------------------------------------------------

ic = 0
DO k=1, sz1
IF (pseudobs_ICRF(k,3) /= 0.d0) THEN
ic = ic +1
pseudobs_ICRF2(ic,1:8) = pseudobs_ICRF(k,1:8)
END IF
END DO


it = 0
DO k=1, sz3
IF (pseudobs_ITRF(k,3) /= 0.d0) THEN
it = it +1
pseudobs_ITRF2(it,1:8) = pseudobs_ITRF(k,1:8)
END IF
END DO
pseudobs_ICRF = 0.d0
pseudobs_ITRF = 0.d0

pseudobs_ICRF = pseudobs_ICRF2
pseudobs_ITRF = pseudobs_ITRF2

END
