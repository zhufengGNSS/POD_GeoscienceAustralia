SUBROUTINE scan0orbext 


! ----------------------------------------------------------------------
! SUBROUTINE: scan0orbext.f03
! ----------------------------------------------------------------------
! Purpose:
!  To skip the bad orbit with zero value for the orbit comparison 
! ----------------------------------------------------------------------
! Input arguments:
! 
! - orbext_ICRF: 	Orbit array (Nx8) in inertial frame (ICRF) including the state vector per epoch
! 					  Collumns desciption per epoch:
!               	- Modified Julian Day number (including the fraction of the day) 
!					- Seconds since 00h 
!					- Position vector (m)
!					- Velocity vector (m/sec)
! - orbext_ITRF:	Orbit array (Nx8) in terrestrial frame (ITRF)
! 					  Collumns desciption per epoch: similar to the orbext_ICRF
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
      INTEGER (KIND = prec_int8) :: sz1, sz2, sz3, sz4, sz5, sz6 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, read_i, k
      INTEGER (KIND = prec_int8) :: ic, it

sz1 = size(orbext_ICRF, DIM = 1)
sz2 = size(orbext_ICRF, DIM = 2)

sz3 = size(orbext_ITRF, DIM = 1)
sz4 = size(orbext_ITRF, DIM = 2)


sz5 = size(orbext_kepler, DIM = 1)
sz6 = size(orbext_kepler, DIM = 2)

!print*,'sz1, sz2 =', sz1, sz2, 'sz3, sz4 =',  sz3, sz4, 'sz5, sz6 =', sz5, sz6
 

orbext_ICRF2  = orbext_ICRF
orbext_ITRF2  = orbext_ITRF

orbext_ICRF2  = 0.d0
orbext_ITRF2  = 0.d0


ic = 0
DO k=1, sz1
IF (orbext_ICRF(k,3) /= 0.d0) THEN
ic = ic +1
orbext_ICRF2(ic,1:8) = orbext_ICRF(k,1:8)
END IF
END DO

it = 0
DO k=1, sz3
IF (orbext_ITRF(k,3) /= 0.d0) THEN
it = it +1
orbext_ITRF2(it,1:8) = orbext_ITRF(k,1:8)

END IF
END DO
orbext_ICRF = 0.d0
orbext_ITRF = 0.d0

orbext_ICRF = orbext_ICRF2
orbext_ITRF = orbext_ITRF2

END
