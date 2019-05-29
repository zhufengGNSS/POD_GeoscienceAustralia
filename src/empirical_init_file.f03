SUBROUTINE empirical_init_file (Niteration, Bias_0, CPR_CS_0)


! ----------------------------------------------------------------------
! SUBROUTINE: empirical_init_file.f03
! ----------------------------------------------------------------------
! Purpose:
!  Orbit empirical parameters initialization 
! ----------------------------------------------------------------------
! Input arguments:
! - Niteration: Orbit parameter estimator's number of iteration 
! - Bias_0: 	Bias parameters acceleration values 
! - CPR_CS_0: 	Cycle per revolution parameters (C and S coefficients) 
!
! Output arguments:
!
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, CRC-SI
! Created:	28 May 2019
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      !USE mdl_param
      USE m_write_prmfile_init
      IMPLICIT NONE
	  
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      INTEGER (KIND = prec_int8), INTENT(IN) :: Niteration 
	  REAL (KIND = prec_q), INTENT(IN) :: Bias_0(3)
	  REAL (KIND = prec_q), INTENT(IN) :: CPR_CS_0(3,2)
! ----------------------------------------------------------------------
! OUT

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------  
      INTEGER (KIND = prec_int8) :: i 
      CHARACTER (LEN=100) :: fname, fname1, fname2				
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
      REAL (KIND = prec_d) :: apriori_3(3), apriori_2(2) 				
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
      CHARACTER (LEN=100), DIMENSION(:), ALLOCATABLE :: param_id_array				
      CHARACTER (LEN=500), DIMENSION(:), ALLOCATABLE :: param_value_array		
      INTEGER (KIND = prec_int8) :: Nparam	  
! ----------------------------------------------------------------------  



! ----------------------------------------------------------------------
! Initial conditions of empirical parameters files
! ----------------------------------------------------------------------
!Bias_accel_aposteriori = Bias_0
!CPR_CS_aposteriori = CPR_CS_0

i = Niteration
write (fname_id, *) i

! ----------------------------------------------------------------------
! Empirical parameters
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Bias 
! ----------------------------------------------------------------------
!fname = 'emp_bias.in'
!param_id = 'bias_rtn'
!apriori_3 = Bias_0
!write (param_value, *) apriori_3
!Call write_prmfile (fname, fname_id, param_id, param_value)

!fname = 'emp_bias.in'
!Nparam = 1
!ALLOCATE (param_id_array(Nparam), STAT = AllocateStatus)
!ALLOCATE (param_value_array(Nparam), STAT = AllocateStatus)
!param_id_array(1) = 'bias_rtn'
!!write (param_value_array(1), *), Bias_0(1), Bias_0(2), Bias_0(3)
!write (param_value_array(1), *), Bias_0
!CALL write_prmfile_init (fname, param_id_array, param_value_array)
!DEALLOCATE (param_id_array, STAT = DeAllocateStatus)
!DEALLOCATE (param_value_array, STAT = DeAllocateStatus)

fname = 'emp_bias.in'
param_id = 'bias_rtn'
apriori_3 = Bias_0
write (param_value, *) apriori_3
CALL write_prmfile_init0 (fname, param_id, param_value)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! CPR terms
! ----------------------------------------------------------------------
!fname = 'emp_cpr.in'
!param_id = 'cpr_radial'
!apriori_2 = (/ CPR_CS_0(1,1) , CPR_CS_0(1,2) /)
!write (param_value, *) apriori_2
!Call write_prmfile (fname, fname_id, param_id, param_value)

!param_id = 'cpr_along'
!apriori_2 = (/ CPR_CS_0(2,1) , CPR_CS_0(2,2) /)
!write (param_value, *) apriori_2
!Call write_prmfile (fname, fname_id, param_id, param_value)

!param_id = 'cpr_cross'
!apriori_2 = (/ CPR_CS_0(3,1) , CPR_CS_0(3,2) /)
!write (param_value, *) apriori_2
!Call write_prmfile (fname, fname_id, param_id, param_value)
! ----------------------------------------------------------------------

fname = 'emp_cpr.in'
Nparam = 3
ALLOCATE (param_id_array(Nparam), STAT = AllocateStatus)
ALLOCATE (param_value_array(Nparam), STAT = AllocateStatus)

param_id_array(1) = 'cpr_radial'
param_id_array(2) = 'cpr_along'
param_id_array(3) = 'cpr_cross'
write (param_value_array(1), *), CPR_CS_0(1,1), CPR_CS_0(1,2)
write (param_value_array(2), *), CPR_CS_0(2,1), CPR_CS_0(2,2)
write (param_value_array(3), *), CPR_CS_0(3,1), CPR_CS_0(3,2)
!print *,"param_id_array", param_id_array
!print *,"param_value_array", param_value_array

CALL write_prmfile_init (fname, param_id_array, param_value_array)
DEALLOCATE (param_id_array, STAT = DeAllocateStatus)
DEALLOCATE (param_value_array, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------


END SUBROUTINE

