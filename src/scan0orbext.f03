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
      INTEGER (KIND = prec_int8) :: ic, it, j
      INTEGER (KIND = prec_int8) :: INEW1, INEW2
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE ::orbext_ITRF2,orbext_ICRF2
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE ::orbext_ITRF3,orbext_ICRF3


sz1 = size(orbext_ICRF, DIM = 1)
sz2 = size(orbext_ICRF, DIM = 2)

sz3 = size(orbext_ITRF, DIM = 1)
sz4 = size(orbext_ITRF, DIM = 2)


sz5 = size(orbext_kepler, DIM = 1)
sz6 = size(orbext_kepler, DIM = 2)

!print*,'sz1, sz2 =', sz1, sz2, 'sz3, sz4 =',  sz3, sz4, 'sz5, sz6 =', sz5, sz6
 
ALLOCATE (orbext_ICRF2(sz1,sz2), STAT = AllocateStatus)
ALLOCATE (orbext_ITRF2(sz3,sz4), STAT = AllocateStatus)

!orbext_ICRF2  = orbext_ICRF
!orbext_ITRF2  = orbext_ITRF

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

ALLOCATE (orbext_ICRF3(ic,sz2), STAT = AllocateStatus)
ALLOCATE (orbext_ITRF3(it,sz3), STAT = AllocateStatus)

INEW1 = 0
DO k=1, sz1
   DO j=1, ic
   IF (orbext_ICRF2(j,3) == orbext_ICRF(k,3) .AND. orbext_ICRF2(j,4) == orbext_ICRF(k,4))THEN

   INEW1 = INEW1 + 1
   orbext_ICRF3(INEW1,1:8)= orbext_ICRF(k,1:8)

   END IF
   END DO
END DO


INEW2 = 0
DO k=1, sz3
   DO j=1, it
   IF (orbext_ITRF2(j,3) == orbext_ITRF(k,3) .AND. orbext_ITRF2(j,4) == orbext_ITRF(k,4))THEN

   INEW2 = INEW2 + 1
   orbext_ITRF3(INEW2,1:8)= orbext_ITRF(k,1:8)

   END IF
   END DO
END DO

DEALLOCATE (orbext_ICRF, STAT = DeAllocateStatus)
DEALLOCATE (orbext_ITRF, STAT = DeAllocateStatus)

ALLOCATE (orbext_ICRF(ic,sz2), STAT = AllocateStatus)
ALLOCATE (orbext_ITRF(it,sz3), STAT = AllocateStatus)
print*,'Number of Epochs or Positions used for orbit comparison =', ic


orbext_ICRF = 0.d0
orbext_ITRF = 0.d0

orbext_ICRF = orbext_ICRF3
orbext_ITRF = orbext_ITRF3

END
