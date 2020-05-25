MODULE m_ecom_init

contains

SUBROUTINE ecom_init (ECOM_0_coef)


! ----------------------------------------------------------------------
! SUBROUTINE: ecom_init
! ----------------------------------------------------------------------
! Purpose: Initialize the ECOM model 
! ----------------------------------------------------------------------
! Input arguments:
! - ECOM_0_coef: 	Array for the initialization of the ECOM model  
! ----------------------------------------------------------------------
! Author :	Tzupang Tseng,  Geoscience Australia
!
! Created:	28-05-2019 
!
! Changes:      03-12-2019 Tzupang Tseng: a option for the simple box wing  
!                                         is created
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_config
      USE mdl_param
      IMPLICIT NONE
	  
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
       REAL (KIND = prec_q), DIMENSION(:), INTENT(OUT), ALLOCATABLE :: ECOM_0_coef
! ----------------------------------------------------------------------
!      INTEGER (KIND = prec_int2) :: Niteration,srp_i 
      INTEGER (KIND = prec_int2) :: srp_i 
      INTEGER (KIND = prec_int8) :: i, ii, PD_Param_ID 
      INTEGER (KIND = prec_int2) :: AllocateStatus
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: fname 				
      CHARACTER (LEN=50) :: fname_id				
      CHARACTER (LEN=100) :: param_id				
      CHARACTER (LEN=500) :: param_value				
! ----------------------------------------------------------------------
! Initial conditions for solar radiation pressure
! ----------------------------------------------------------------------
!ALLOCATE (ECOM_0_coef(NPARAM_glb))
ALLOCATE (ECOM_0_coef(ECOMNUM))
srp_i = ECOM_param_glb
!i = Niteration
!write (fname_id, *) i
ECOM_0_coef = 0.d0
IF (IC_MODE_cfg == 2) THEN
!	ECOM_0_coef = IC_sat_glb (9:8+NPARAM_glb)
        ECOM_0_coef = IC_sat_glb (9:8+ECOMNUM)
END IF
!print*,'ECOM_0_coef=',ECOM_0_coef

IF (srp_i == 1) THEN
fname = 'ECOM1_srp.in'
param_id = 'ECOM1'
write (param_value, *) ECOM_0_coef
!Call write_prmfile (fname, fname_id, param_id, param_value)
CALL write_prmfile_init0 (fname, param_id, param_value)
END IF


IF (srp_i == 2) THEN
fname = 'ECOM2_srp.in'
param_id = 'ECOM2'
write (param_value, *) ECOM_0_coef
!Call write_prmfile (fname, fname_id, param_id, param_value)
CALL write_prmfile_init0 (fname, param_id, param_value)
END IF

IF (srp_i == 3) THEN
fname = 'SBOXW_srp.in'
param_id = 'SBOXW'
write (param_value, *) ECOM_0_coef
!Call write_prmfile (fname, fname_id, param_id, param_value)
CALL write_prmfile_init0 (fname, param_id, param_value)
END IF


END SUBROUTINE

END MODULE
